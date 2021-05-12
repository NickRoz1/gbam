use super::{RawRecord, FIELDS_NUM, GBAM_MAGIC};
use crate::meta::RowGroupMeta;
use crate::rowgroup::{CompressionDaemon, RowGroup};
use byteorder::{LittleEndian, WriteBytesExt};
// use flate2::write::ZlibEncoder;
use std::io;
use std::io::Write;
use std::sync::{Arc, Mutex};
// use zstd::stream::Encoder;

/// A GBAM writer.
pub struct Writer {
    cur_block_offset: u64,
    cur_rowgroup: RowGroup,
    // Option utilized to control the drop order.
    compression_daemon: Option<CompressionDaemon>,
    writer_thread: Arc<Mutex<WriterThread>>,
    rowgroups_info: Option<Arc<Mutex<Vec<RowGroupMeta>>>>,
}

/// To leverage column oriented storage the data (iterator) can be requested as
/// special tuple and only the necessary columns will be traversed.
///
/// The rowgroups are of variable size. If variable sized field overflow maximum
/// rowgroup size, it placed in a queue for later processing, and rowgroup
/// without this record is processed.
///
/// Example:
///
///
impl Writer {
    /// Writes ID data to passed writers, wraps them into encoder and creates
    /// Writer object using them.
    pub fn new<T: Write + Send + 'static>(
        mut inner: T,
        mut thread_count: usize,
    ) -> io::Result<Self> {
        // Probably doesnt make sense if more than FIELDS NUM? Rowgroup is a bottleneck currently.
        thread_count = std::cmp::min(num_cpus::get(), thread_count);
        inner.write(&GBAM_MAGIC)?;

        let (buffer_sender, buffer_receiver): (flume::Sender<Vec<u8>>, flume::Receiver<Vec<u8>>) =
            flume::unbounded();

        // There should be at least one buffer per column. If the the last
        // column chunk which will come out out of compressor worker will be the
        // first one out of rowgroup and there are less buffers than column
        // chunks, the sorter thread will keep all the available buffers,
        // blocking daemon thread from loading this first task to unload the
        // sorter, thus introducing a deadlock.
        for _ in 0..FIELDS_NUM {
            // Prefill with buffers
            buffer_sender.send(Vec::<u8>::new()).unwrap();
        }

        let writer_thread = Arc::new(Mutex::new(WriterThread::new(inner, buffer_sender)));

        let rowgroup_meta_vector = Arc::new(Mutex::new(Vec::new()));

        let compression_thread = CompressionDaemon::new(
            rowgroup_meta_vector.clone(),
            buffer_receiver,
            writer_thread.clone(),
            thread_count,
        );

        Ok(Self {
            cur_block_offset: 0,
            cur_rowgroup: RowGroup::default(),
            compression_daemon: Some(compression_thread),
            writer_thread: writer_thread,
            rowgroups_info: Some(rowgroup_meta_vector),
        })
    }

    /// Splits the record field by field and writes them to separate files.
    /// Works correctly only when full record is provided.
    pub fn write_record(&mut self, rec: &RawRecord) -> io::Result<()> {
        if self.cur_rowgroup.payload_size_of() + rec.len()
            > self.cur_rowgroup.max_uncompr_size as usize
        {
            self.flush();
        }
        self.cur_rowgroup.write_record(rec)?;
        Ok(())
    }

    /// Flush the buffer.
    pub fn flush(&mut self) {
        if self.cur_rowgroup.is_empty() {
            return;
        }
        // Receive already dumped rowgroup buffer and swap it with the
        // current one to avoid allocation.
        let daemon = self.compression_daemon.as_mut().unwrap();
        let mut used_buffer = daemon.recv_buf().unwrap();
        std::mem::swap(&mut used_buffer, &mut self.cur_rowgroup);
        // Send rowgroup to compression daemon from wherein it will be sent
        // to IO thread and written into the file
        daemon.compress(used_buffer).unwrap();
        self.cur_rowgroup.clear();
    }

    /// Writes meta into the file. Should only be called before Writer destruction.
    fn write_meta(&mut self) -> io::Result<usize> {
        let file_meta = Arc::try_unwrap(self.rowgroups_info.take().unwrap())
            .unwrap()
            .into_inner()
            .unwrap();
        let mut output_writer = self.writer_thread.lock().unwrap();
        let mut meta_size: usize = 0;
        for rowgroup_meta in file_meta.into_iter() {
            let serialized: Vec<u8> = rowgroup_meta.into();
            meta_size += serialized.len();
            output_writer.send(serialized).unwrap();
        }
        Ok(meta_size)
    }
}

impl Drop for Writer {
    fn drop(&mut self) {
        self.flush();
        // Drop the daemon
        self.compression_daemon = None;
        let meta_size = self.write_meta().unwrap();
        ///// Write additional data + magic string
        let mut output_writer = self.writer_thread.lock().unwrap();
        let buf: Vec<u8> = vec![0; 8 + GBAM_MAGIC.len()];
        let mut writer = std::io::Cursor::new(buf);
        writer.write_u64::<LittleEndian>(meta_size as u64).unwrap();
        writer.write(&GBAM_MAGIC[..]).unwrap();
        output_writer.send(writer.into_inner()).unwrap();
    }
}

/// This struct manages background thread writer.
pub struct WriterThread {
    thread_handle: Option<std::thread::JoinHandle<()>>,
    buf_pipe_t: Option<flume::Sender<Vec<u8>>>,
}

impl WriterThread {
    /// The writer thread is used to avoid using main thread time on IO
    /// operation.
    /// # Input parameters
    /// inner : inner writer, which can be retrieved upon request.
    /// buffer_sender: used to send buffers with already written data back. It is
    /// utilized (instead of incapsulating it inside WriterThread object) to
    /// avoid using mutex in tight loops in compression daemon thread (channel
    /// is preferred).
    pub fn new<T: Write + Send + 'static>(
        mut inner: T,
        buffer_sender: flume::Sender<Vec<u8>>,
    ) -> Self {
        let (control_channel_t, control_channel_r): (
            flume::Sender<Vec<u8>>,
            flume::Receiver<Vec<u8>>,
        ) = flume::unbounded();

        let thread_handle = std::thread::spawn(move || {
            while let Ok(buf) = control_channel_r.recv() {
                if inner.write(&buf[..]).expect("Writer failed.") != buf.len() {
                    panic!("Writing to file failed.");
                }
                // When writing meta the buffers are useless, so drop them (plus there is no receiver by this time).
                if !buffer_sender.is_disconnected() {
                    buffer_sender.send(buf).unwrap();
                }
            }
        });

        WriterThread {
            thread_handle: Some(thread_handle),
            buf_pipe_t: Some(control_channel_t),
        }
    }

    /// Sends data to writer thread to write into file
    pub fn send(&mut self, buf: Vec<u8>) -> Result<(), flume::SendError<Vec<u8>>> {
        self.buf_pipe_t.as_ref().unwrap().send(buf)
    }
}

impl Drop for WriterThread {
    fn drop(&mut self) {
        // Drop the channel
        self.buf_pipe_t.take();
        self.thread_handle
            .take()
            .unwrap()
            .join()
            .expect("Writer thread failed.");
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;
    #[test]
    fn test_writer_thread() {
        static mut BUF: [u8; 5] = [0, 0, 0, 0, 0];
        let data: Vec<u8> = vec![0, 1, 2, 3, 4];
        let (buffer_sender, buffer_receiver): (flume::Sender<Vec<u8>>, flume::Receiver<Vec<u8>>) =
            flume::unbounded();

        unsafe {
            let out = Cursor::new(&mut BUF[..]);

            let mut writer_thread = WriterThread::new(out, buffer_sender);
            writer_thread.send(data.clone()).unwrap();
            assert_eq!(data, buffer_receiver.recv().unwrap());
            drop(writer_thread);
            assert_eq!(&data[..], BUF);
        }
    }

    #[test]
    fn test_writer() {
        let records = vec![RawRecord::default(); 4];
        static mut BUF: [u8; 10000] = [0; 10000];

        unsafe {
            let out = Cursor::new(&mut BUF[..]);

            let mut writer = Writer::new(out, 4).unwrap();
            writer.write_record(&records[0]);
            writer.flush();
            // for record in records.iter() {
            //     writer.write_record(record);
            // }
            // writer.flush();

            drop(writer);

            for i in (0..10000).rev() {
                if BUF[i] != 0 {
                    println!("THE OCCUPIED LENGTH = {:?}", i);
                    break;
                }
            }
        }
    }
}
