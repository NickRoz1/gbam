// This structure splits records and fills column buffers, which later get
// written into the output stream.

use super::{mega_byte_size, ColChunkMeta, RowGroupMeta};
use super::{Compression, Fields, RawRecord, FIELDS_NUM};
use crate::writer::WriterThread;
use byteorder::{LittleEndian, ReadBytesExt, WriteBytesExt};
use flate2::write::ZlibEncoder;
use flume;
use std::convert::TryInto;
use std::io::Write;
use std::sync::{Arc, Mutex};
use zstd::Encoder;

pub struct RowGroup {
    offsets: Vec<u64>,
    // Compression type, number of values, bytes of column
    columns: Vec<(Compression, usize, Vec<u8>)>,
    // If the whole records are queried, the decompressed columns of rowgroup will be smaller than this size.
    pub max_uncompr_size: u64,
}

impl RowGroup {
    pub fn write_record(&mut self, rec: &RawRecord) -> std::io::Result<()> {
        for field in Fields::iterator() {
            match field {
                Fields::LName | Fields::SequenceLength | Fields::NCigar | Fields::RawTagsLen => {
                    let offset = &mut self.offsets[*field as usize];
                    match field {
                        // Columns, containing offsets (total) to the variable sized fields beginnings.
                        // So the offset to third Sequence equals to 0 + Seq_1.len() + Seq_2.len();
                        // So the third Sequence goes after [0..Seq_1.len()] and [Seq_1.len()..Seq_2.len()]
                        // till the end of the column.
                        Fields::LName => self.columns[*field as usize].2.write_u8(*offset as u8)?,
                        Fields::SequenceLength => self.columns[*field as usize]
                            .2
                            .write_u32::<LittleEndian>(*offset as u32)?,
                        Fields::NCigar => self.columns[*field as usize]
                            .2
                            .write_u16::<LittleEndian>(*offset as u16)?,
                        Fields::RawTagsLen => self.columns[*field as usize]
                            .2
                            .write_u32::<LittleEndian>(*offset as u32)?,
                        _ => panic!("This field is not supported: {} \n", *field as usize),
                    }
                    *offset += rec.get_len_val(field) as u64;
                }
                _ => {
                    // Variable sized fields
                    self.columns[*field as usize]
                        .2
                        .write(rec.get_bytes(field))?;
                }
            }
            self.columns[*field as usize].1 += 1;
        }

        Ok(())
    }

    /// Use this function to customize compression codecs for columns of this rowgroup.
    pub fn compression_policy(&mut self, compression_types: Vec<Compression>) {
        for (col, compr_type) in self.columns.iter_mut().zip(compression_types) {
            col.0 = compr_type;
        }
    }

    /// Clears payload containing fields
    pub fn clear(&mut self) {
        for offset in self.offsets.iter_mut() {
            *offset = 0;
        }
        for col in self.columns.iter_mut() {
            col.1 = 0;
            col.2.clear();
        }
    }

    pub fn is_empty(&self) -> bool {
        self.columns.iter().all(|(_, num_vals, _)| *num_vals == 0)
    }

    // Get size of actual payload
    pub fn payload_size_of(&self) -> usize {
        let mut total_size: usize = 0;
        for (_, _, col_payload) in self.columns.iter() {
            total_size += col_payload.len();
        }
        total_size
    }

    /// Get byte size of this rowgroup (compression and number of values does
    /// not go to the disk so even if the rowgroup is empty it's columns are
    /// not!).
    pub fn size_of(&self) -> usize {
        let mut total_size: usize = 0;
        for col in self.columns.iter() {
            total_size += std::mem::size_of_val(&*self.offsets);
        }
        total_size
    }
}

impl Default for RowGroup {
    fn default() -> Self {
        RowGroup {
            offsets: vec![0; FIELDS_NUM],
            columns: vec![(Compression::default(), 0, Vec::new()); FIELDS_NUM],
            max_uncompr_size: 16 * mega_byte_size as u64,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_rowgroup() -> () {
        let mut RowGroup = RowGroup::default();
        let mut record = RawRecord::default();
        RowGroup.write_record(&record).unwrap();
        let mut assemled_record: Vec<u8> = Vec::<u8>::new();
        for mut col in RowGroup.columns.iter_mut() {
            assemled_record.append(&mut col.2);
        }
        // Nullify name length since the RowGroup counts total offsets in the
        // variable sized fields lengths columns. Offset to first value = 0.
        record.0[8] = 0;
        // Omit RawTagsLen (u32) since it's not included in the original record.
        assert_eq!(record.0[..], assemled_record[0..assemled_record.len() - 4]);
    }
}

// In main thread while let on the reader pipe, each blob from reader gets
// parsed and grouped into rowgroup, then rowgroup is sent ot daemon thread,
// daemon thread dispenses columns chunks between worker threads, after
// compression it receives the compressed chunks, "places" them in order, dumps
// it into writer and creates rowgroup meta for it.
//
// The rowgroup meta vector (containing all rowgroups meta objects for this
// input) is shared with daemon thread via ARC MUTEX. After joining daemon thread
// and writter thread, the rowgroup meta vector will contain meta information
// about written data. GBAM file has been created.

/// Column chunk compression task struct
struct CompressionTask {
    buf_num: Option<usize>,
    compr_type: Compression,
    val_num: usize,
    payload: Vec<u8>,
    compressed: Vec<u8>,
}

/// This struct responsible for compression in separate thread.
pub struct CompressionDaemon {
    // Reason for option: https://users.rust-lang.org/t/spawn-threads-and-join-in-destructor/1613/2
    daemon_thread_handle: Option<std::thread::JoinHandle<()>>,
    workers_handles: Vec<Option<std::thread::JoinHandle<()>>>,
    writer_thread: Arc<Mutex<WriterThread>>,
    ordering_thread_handle: Option<std::thread::JoinHandle<()>>,
    // Used to send rowgroups to compress to daemon thread. Option utilized to
    // be able to manually control drop order in the destructor of Compression
    // Daemon.
    compression_channel: Option<flume::Sender<RowGroup>>,
    processed_rowgroups_channel: Option<flume::Receiver<RowGroup>>,
}

/// Output buffer. Usize used to order compressed chunks.
#[derive(PartialEq, Eq)]
struct OutputBuf(usize, Vec<u8>);

impl Ord for OutputBuf {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // Flipped comparison call to transform max heap into min heap
        other.0.cmp(&self.0)
    }
}
impl PartialOrd for OutputBuf {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl CompressionDaemon {
    /// The compression daemon performs compression tasks in separate threads.
    /// # Input parameters
    ///
    /// * **rowgroups_meta**: vector containing meta information about rowgroups.
    /// * **buffer_receiver**: receives drained buffers from output performing thread.    
    /// * **thread_count**: amount of threads to use for compression, but no less than 2. (Affects memory usage)  
    pub fn new(
        rowgroup_meta_vector: Arc<Mutex<Vec<RowGroupMeta>>>,
        buffer_receiver: flume::Receiver<Vec<u8>>,
        writer_thread: Arc<Mutex<WriterThread>>,
        thread_count: usize,
    ) -> Self {
        assert!(thread_count >= 3);

        let mut worker_handles = Vec::<Option<std::thread::JoinHandle<()>>>::new();

        // Channel used to return compressed rowgroups back to supplier thread.
        let (processed_rowgroups_t, processed_rowgroups_r): (
            flume::Sender<RowGroup>,
            flume::Receiver<RowGroup>,
        ) = flume::unbounded();

        // Channel used to supply daemon thread with rowgroups to compress.
        let (rowgroup_pipe_t, rowgroup_pipe_r): (
            flume::Sender<RowGroup>,
            flume::Receiver<RowGroup>,
        ) = flume::unbounded();

        // This channel used to supply tasks to worker threads.
        let (tasks_sink_t, tasks_stream_r): (
            flume::Sender<CompressionTask>,
            flume::Receiver<CompressionTask>,
        ) = flume::unbounded();

        // Used to return finished tasks back to supplying thread.
        let (finished_tasks_t, finished_tasks_stream_r): (
            flume::Sender<CompressionTask>,
            flume::Receiver<CompressionTask>,
        ) = flume::unbounded();

        // Used to send compressed data to ordering thread which ensures that column chunks written in order.
        let (sorter_t, sorter_buffers_r): (flume::Sender<OutputBuf>, flume::Receiver<OutputBuf>) =
            flume::unbounded();

        for _ in 0..(thread_count - 2) {
            let tasks_stream_r_clone = tasks_stream_r.clone();
            let finished_tasks_t_clone = finished_tasks_t.clone();
            let compr_task = CompressionTask {
                buf_num: None,
                compr_type: Compression::ZSTD(0),
                val_num: 0,
                payload: Vec::new(),
                compressed: Vec::new(),
            };
            // Prefill the queue with "completed" tasks to initiate processing
            finished_tasks_t.send(compr_task).unwrap();
            worker_handles.push(Some(std::thread::spawn(move || {
                while let Ok(mut task) = tasks_stream_r_clone.recv() {
                    // println!("IN TASK");
                    task.compressed.clear();
                    let compressed_bytes = match task.compr_type {
                        Compression::FLATE2(level) => {
                            let mut e = ZlibEncoder::new(
                                task.compressed,
                                flate2::Compression::new((level as u32).try_into().unwrap()),
                            );
                            e.write(&task.payload);
                            e.finish()
                        }
                        Compression::ZSTD(level) => {
                            let mut e = zstd::stream::Encoder::new(
                                task.compressed,
                                (level as u32).try_into().unwrap(),
                            )
                            .unwrap();
                            // println!("IN TASK: {:?} / {:?}", task.payload, e.write(&task.payload));
                            e.finish()
                        }
                    };
                    task.compressed = compressed_bytes.unwrap();
                    println!("Leaving TASK");
                    println!("PROCESSED TASK {}", task.buf_num.unwrap());
                    finished_tasks_t_clone.send(task).unwrap();
                }
                println!("TASK THREAD HAS BEEN DROPPED");
            })));
        }

        let writer_thread_clone = writer_thread.clone();
        // This thread orders compressed column chunks before dumping them into I/O thread (writer thread).
        let ordering_thread_handle = Some(std::thread::spawn(move || {
            use std::collections::BinaryHeap;
            let mut bin_heap = BinaryHeap::<OutputBuf>::new();
            let mut cur_chunk_num: usize = 0;
            let mut output = writer_thread_clone.lock().unwrap();
            while let Ok(output_buf) = sorter_buffers_r.recv() {
                bin_heap.push(output_buf);
                println!("IN ORDERING THREAD");
                while !bin_heap.is_empty() && bin_heap.peek().unwrap().0 == cur_chunk_num {
                    output.send(bin_heap.pop().unwrap().1).unwrap();
                    cur_chunk_num += 1;
                }
                if !bin_heap.is_empty() {
                    println!(
                        "AWAITING IN HEAP {} CUR TOP {}",
                        bin_heap.len(),
                        bin_heap.peek().unwrap().0
                    );
                }
            }
            while !bin_heap.is_empty() {
                assert_eq!(bin_heap.peek().unwrap().0, cur_chunk_num);
                output.send(bin_heap.pop().unwrap().1).unwrap();
                cur_chunk_num += 1;
            }
        }));

        let processed_data_t = processed_rowgroups_t.clone();
        let thread_handle = Some(std::thread::spawn(move || {
            // Offset into the current writer.
            let mut writer_offset: u64 = 0;

            let mut cur_buf = 0;
            let mut processed_bufs_num = 0;
            // Get rowgroup from supplying thread.
            while let Ok(mut rowgroup) = rowgroup_pipe_r.recv() {
                println!("CUR ROWGROUP LEN - {:?}", rowgroup.columns.len());
                let mut rowgroup_meta = RowGroupMeta {
                    offset: writer_offset,
                    cols: vec![ColChunkMeta::default(); FIELDS_NUM],
                };

                let mut relative_offset: u64 = 0;

                for col_to_fill_idx in 0..FIELDS_NUM {
                    let mut buf = finished_tasks_stream_r.recv().unwrap();
                    // If this buffer wasn't utilized yet (it is prefilled).
                    if buf.buf_num.is_some() {
                        processed_bufs_num += 1;
                        println!("Got the task number {}", buf.buf_num.unwrap());
                        // The compressed buffers are written in the same order as
                        // they came in, and the reader of BAM files orders the
                        // records too.
                        let col_num = buf.buf_num.unwrap() % FIELDS_NUM;
                        rowgroup_meta.cols[col_num] = ColChunkMeta {
                            uncompr_size: buf.payload.len() as u64,
                            val_num: buf.val_num as u64,
                            compr_size: buf.compressed.len() as u64,
                            offset: writer_offset + relative_offset,
                            compressor: buf.compr_type,
                        };
                        // Incorrect! Doesnt account for order! Task may come out of order.
                        relative_offset += buf.compressed.len() as u64;
                        {
                            println!("DEADLOCK IS");
                            let mut output_buf = buffer_receiver.recv().unwrap();
                            println!("DEADLOCK IS HERE");
                            std::mem::swap(&mut output_buf, &mut buf.compressed);
                            sorter_t
                                .send(OutputBuf(buf.buf_num.unwrap(), output_buf))
                                .unwrap();
                            println!("GOT BACK");
                        }
                    }
                    let (compr, val_num, col) = &mut rowgroup.columns[col_to_fill_idx];
                    // Load the buffer with new data
                    {
                        std::mem::swap(col, &mut buf.payload);
                        buf.compr_type = compr.clone();
                        buf.buf_num = Some(cur_buf);
                        buf.val_num = *val_num;
                        // Assumed the size of compressed data won't exceed the
                        // size of uncompressed data by more than 200 bytes.
                        buf.compressed.clear();
                        cur_buf += 1;
                    }
                    tasks_sink_t.send(buf).unwrap();
                }
                println!("FINISHED");
                writer_offset += rowgroup_meta
                    .cols
                    .iter()
                    .fold(0, |acc, col| acc + col.compr_size);
                rowgroup_meta_vector.lock().unwrap().push(rowgroup_meta);
                processed_data_t.send(rowgroup).unwrap();
            }

            let mut rowgroup_meta = RowGroupMeta {
                offset: writer_offset,
                cols: vec![ColChunkMeta::default(); FIELDS_NUM],
            };
            println!("HERE");
            let mut relative_offset: u64 = 0;
            // Process data which is still in compression workers queues.
            loop {
                if processed_bufs_num == cur_buf {
                    // Compression has been finished, free the workers.
                    println!("THE DAEMON EXITING");
                    drop(tasks_sink_t);
                    break;
                }
                let mut buf = finished_tasks_stream_r.recv().unwrap();
                // If this buffer wasn't utilized yet (it is prefilled) dot.
                if buf.buf_num.is_some() {
                    processed_bufs_num += 1;
                    println!("TERMINATING : {}", buf.buf_num.unwrap());
                    // The compressed buffers are written in the same order as
                    // they came in, and the reader of BAM files orders the
                    // records too.
                    let col_num = buf.buf_num.unwrap() % FIELDS_NUM;
                    rowgroup_meta.cols[col_num] = ColChunkMeta {
                        uncompr_size: buf.payload.len() as u64,
                        val_num: buf.val_num as u64,
                        compr_size: buf.compressed.len() as u64,
                        offset: writer_offset + relative_offset,
                        compressor: buf.compr_type,
                    };
                    relative_offset += buf.compressed.len() as u64;
                    {
                        let mut output_buf = buffer_receiver.recv().unwrap();
                        std::mem::swap(&mut output_buf, &mut buf.compressed);
                        sorter_t
                            .send(OutputBuf(buf.buf_num.unwrap(), output_buf))
                            .unwrap();
                        println!("GOT BACK: {}", buf.buf_num.unwrap());
                    }
                }
            }
            rowgroup_meta_vector.lock().unwrap().push(rowgroup_meta);
        }));

        // Prefill the rowgroup buffers queue
        for _ in 0..10 {
            processed_rowgroups_t.send(RowGroup::default()).unwrap();
        }
        CompressionDaemon {
            daemon_thread_handle: thread_handle,
            ordering_thread_handle: ordering_thread_handle,
            workers_handles: worker_handles,
            writer_thread: writer_thread,
            compression_channel: Some(rowgroup_pipe_t),
            processed_rowgroups_channel: Some(processed_rowgroups_r),
        }
    }
    /// Receive used buffer from compression daemon.
    ///
    /// This is a blocking operation. If the buffer is not yet processed, the
    /// thread from where this function where called will be locked until the
    /// buffer is processed.
    pub fn recv_buf(&mut self) -> Result<RowGroup, flume::RecvError> {
        self.processed_rowgroups_channel.as_ref().unwrap().recv()
    }

    /// Send buffer to compression daemon.
    pub fn compress(&mut self, buf: RowGroup) -> Result<(), flume::SendError<RowGroup>> {
        self.compression_channel.as_ref().unwrap().send(buf)
    }
}

impl Drop for CompressionDaemon {
    fn drop(&mut self) {
        // Drop the sender (will lead to daemon thread termination)
        self.compression_channel = None;
        println!("DESTRUCTOR CALLED");
        self.daemon_thread_handle
            .take()
            .unwrap()
            .join()
            .expect("Daemon thread failed.");
        println!("DESTRUCTOR CALLED");
        for handle in self.workers_handles.iter_mut() {
            handle
                .take()
                .unwrap()
                .join()
                .expect("Worker thread failed.");
        }
        self.ordering_thread_handle
            .take()
            .unwrap()
            .join()
            .expect("Ordering thread failed.");
        assert!(!self.writer_thread.is_poisoned());
    }
}
