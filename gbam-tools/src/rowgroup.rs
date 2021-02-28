// This structure splits records and fills column buffers, which later get
// written into the output stream.

use super::{ColChunkMeta, RowGroupMeta};
use super::{Compression, Fields, RawRecord, FIELDS_NUM};
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
    max_uncompr_size: u64,
    cur_size: u64,
    out_buffer: Vec<u8>,
}

impl RowGroup {
    fn write_record(&mut self, rec: &RawRecord) -> std::io::Result<()> {
        if self.cur_size + rec.len() as u64 > self.max_uncompr_size {
            self.flush();
            self.write_record(rec)?;
        }

        for field in Fields::iterator() {
            match field {
                Fields::LName | Fields::SequenceLength | Fields::NCigar | Fields::RawTagsLen => {
                    let offset = &mut self.offsets[*field as usize];
                    match field {
                        // Fixed sized fields
                        Fields::LName => self.columns[*field as usize]
                            .write_u32::<LittleEndian>(*offset as u32)?,
                        Fields::SequenceLength => self.columns[*field as usize]
                            .write_u32::<LittleEndian>(*offset as u32)?,
                        Fields::NCigar => self.columns[*field as usize]
                            .write_u32::<LittleEndian>(*offset as u32)?,
                        Fields::RawTagsLen => self.columns[*field as usize]
                            .write_u32::<LittleEndian>(*offset as u32)?,
                        _ => panic!("This field is not supported: {} \n", *field as usize),
                    }
                    *offset += rec.get_len_val(field) as u64;
                }
                _ => {
                    // Variable sized fields
                    self.columns[*field as usize].write(rec.get_bytes(field))?;
                }
            }
        }

        Ok(())
    }

    fn flush(&mut self) -> Result<()> {}
}

impl Drop for RowGroup {
    fn drop(&mut self) {
        self.flush();
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
    buf_num: usize,
    compr_type: Compression,
    payload: Vec<u8>,
    compressed: Vec<u8>,
}

/// This struct responsible for compression in separate thread.
struct CompressionDaemon<T> {
    daemon_thread_handle: std::thread::JoinHandle<()>,
    workers_handles: Vec<std::thread::JoinHandle<()>>,
    writer_thread: WriterThread<T>,
    ordering_thread_handle: std::thread::JoinHandle<()>,
    // Used to send rowgroups to compress to daemon thread.
    compression_channel: flume::Sender<RowGroup>,
}

/// Output buffer. Usize used to order compressed chunks.
struct OutputBuf(usize, Vec<u8>);

impl Ord for OutputBuf {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // Flipped comparison call to transform max heap into min heap
        other.0.cmp(&self.0)
    }
}

impl<T> CompressionDaemon<T>
where
    T: Write + Send + 'static,
{
    /// The compression daemon performs compression tasks in separate threads.
    /// # Input parameters
    /// buffer_pipe: supplies buffers for compression.
    /// buffer_receiver: receives drained buffers from output performing thread.
    /// thread_count: amount of threads to use for compression, but no less than 2. (Affects memory usage)
    pub fn new(
        rowgroup_buffer_pipe: flume::Receiver<RowGroup>,
        rowgroups_meta: Arc<Mutex<Vec<RowGroupMeta>>>,
        buffer_receiver: flume::Receiver<Vec<u8>>,
        writer_thread: WriterThread<T>,
        thread_count: usize,
    ) -> (Self, flume::Receiver<RowGroup>) {
        assert!(thread_count >= 2);

        let mut worker_handles = Vec::<std::thread::JoinHandle<()>>::new();

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

        // Used to return drained buffers back to worker threads.
        let (dumped_buffers_sink, dumped_buffers_stream_r): (
            flume::Sender<Vec<u8>>,
            flume::Receiver<Vec<u8>>,
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

        for _ in 0..(thread_count - 1) {
            let tasks_stream_r_clone = tasks_stream_r.clone();
            worker_handles.push(std::thread::spawn(move || {
                while let Ok(mut task) = tasks_stream_r_clone.recv() {
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
                            e.write(&task.payload);
                            e.finish()
                        }
                    };
                    finished_tasks_t.send(task).unwrap();
                }
            }));
        }

        // This thread orders compressed column chunks before dumping them into I/O thread (writer thread).
        let ordering_thread_handle = std::thread::spawn(move || {
            use std::collections::BinaryHeap;
            let mut bin_heap = BinaryHeap::<OutputBuf>::new();
            let mut cur_chunk_num: usize = 0;
            let writer_thread = writer_thread;
            while let Ok(output_buf) = sorter_buffers_r.recv() {
                bin_heap.push(output_buf);
                while !bin_heap.is_empty() && bin_heap.peek().unwrap().0 == cur_chunk_num {
                    writer_thread
                        .buf_pipe_t
                        .send(bin_heap.pop().unwrap().1)
                        .unwrap();
                    cur_chunk_num += 1;
                }
            }
            while !bin_heap.is_empty() {
                assert_eq!(bin_heap.peek().unwrap().0, cur_chunk_num);
                writer_thread
                    .buf_pipe_t
                    .send(bin_heap.pop().unwrap().1)
                    .unwrap();
                cur_chunk_num += 1;
            }
            writer_thread.thread_handle.join().unwrap();
        });

        let thread_handle = std::thread::spawn(move || {
            let rowgroup_meta_vector = rowgroups_meta.clone();
            // Offset into the current writer.
            let writer_offset: u64 = 0;

            let cur_buf = 0;
            // Get rowgroup from supplying thread.
            while let Ok(rowgroup) = rowgroup_pipe_r.recv() {
                let rowgroup_meta = RowGroupMeta {
                    offset: writer_offset,
                    cols: vec![ColChunkMeta::default(); FIELDS_NUM],
                };

                let relative_offset: u64 = 0;

                for (compr, _, col) in rowgroup.columns.iter_mut() {
                    let buf = finished_tasks_stream_r.recv().unwrap();
                    let col_num = buf.buf_num % FIELDS_NUM;
                    rowgroup_meta.cols[col_num] = ColChunkMeta {
                        uncompr_size: rowgroup.columns[col_num].2.len() as u64,
                        val_num: rowgroup.columns[col_num].1 as u64,
                        compr_size: buf.compressed.len() as u64,
                        offset: writer_offset + relative_offset,
                        compressor: buf.compr_type,
                    };

                    relative_offset += buf.compressed.len() as u64;
                    // Sink the buffer
                    {
                        let output_buf = buffer_receiver.recv().unwrap();
                        std::mem::swap(&mut output_buf, &mut buf.compressed);
                        sorter_t.send(OutputBuf(buf.buf_num, output_buf)).unwrap();
                    }
                    // Load the buffer with new data
                    {
                        std::mem::swap(col, &mut buf.payload);
                        buf.compr_type = *compr;
                        buf.buf_num = cur_buf;
                        cur_buf += 1;
                    }

                    tasks_sink_t.send(buf).unwrap();
                }

                writer_offset += rowgroup_meta
                    .cols
                    .iter()
                    .fold(0, |acc, col| acc + col.compr_size);
                processed_rowgroups_t.send(rowgroup).unwrap();
            }
        });

        (
            CompressionDaemon::<T> {
                daemon_thread_handle: thread_handle,
                ordering_thread_handle: ordering_thread_handle,
                workers_handles: worker_handles,
                writer_thread: writer_thread,
                compression_channel: rowgroup_pipe_t,
            },
            processed_rowgroups_r,
        )
    }
    // Receive used buffer from compression daemon
    pub fn recv() -> Vec<u8> {}

    pub fn compress(buf: Vec<u8>) {}
}

// This struct manages background thread writer.
struct WriterThread<T: Write> {
    thread_handle: std::thread::JoinHandle<()>,
    buf_pipe_t: flume::Sender<Vec<u8>>,
}

impl<T> WriterThread<T>
where
    T: Write + Send + 'static,
{
    /// The writer thread is used to avoid using main thread time on IO operation.
    /// # Input parameters
    /// inner : inner writer, which can be retrieved upon request.
    /// buffer_sender: used to retrieve buffers with already written data back.
    pub fn new(inner: T, buffer_sender: flume::Sender<Vec<u8>>) -> Self {
        let (control_channel_t, control_channel_r): (
            flume::Sender<Vec<u8>>,
            flume::Receiver<Vec<u8>>,
        ) = flume::unbounded();

        let thread_handle = std::thread::spawn(move || {
            while let Ok(buf) = control_channel_r.recv() {
                if inner.write(&buf[..]).expect("Writer failed.") != buf.len() {
                    panic!("Writing to file failed.");
                }
                buffer_sender.send(buf).unwrap();
            }
        });

        WriterThread::<T> {
            thread_handle: thread_handle,
            buf_pipe_t: control_channel_t,
        }
    }
}

impl<T> Drop for WriterThread<T> {
    fn drop(&mut self) {
        self.thread_handle.join();
    }
}
