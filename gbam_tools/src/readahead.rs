// This module preparses GBAM blocks to parallelize decompression
use super::{decompress_block, Codecs, ReadSeekSendStatic};
use flume::{Receiver, Sender};
use rayon::ThreadPool;
use std::collections::VecDeque;
use std::io::{Read, Seek, SeekFrom};
use std::sync::{Arc, Condvar, Mutex};
/// Prefetches and decompresses GBAM blocks
pub(crate) struct Readahead {
    inner: Arc<Mutex<Box<dyn ReadSeekSendStatic>>>,
    pool: ThreadPool,
    read_bufs_queue: Arc<Mutex<VecDeque<Vec<u8>>>>,
    sender: Sender<Vec<u8>>,
    receiver: Receiver<Vec<u8>>,
    // Used to ensure the uncompressed blocks go in request order
    cur_task: usize,
    cvar_task_num: Arc<(Mutex<usize>, Condvar)>,
}

impl Readahead {
    pub fn new(thread_num: usize, reader: Box<dyn ReadSeekSendStatic>) -> Self {
        let (tx, rx) = flume::unbounded();
        let mut read_bufs_queue = VecDeque::with_capacity(thread_num);
        (0..thread_num).for_each(|_| read_bufs_queue.push_back(Vec::new()));
        Self {
            pool: rayon::ThreadPoolBuilder::new()
                .num_threads(thread_num)
                .build()
                .unwrap(),
            inner: Arc::new(Mutex::new(reader)),
            read_bufs_queue: Arc::new(Mutex::new(read_bufs_queue)),
            sender: tx,
            receiver: rx,
            cur_task: 0,
            cvar_task_num: Arc::new((Mutex::new(0), Condvar::new())),
        }
    }

    /// Creates task to parse block at position of certain size
    pub fn create_task(
        &mut self,
        seek_pos: u64,
        block_size: usize,
        mut buf: Vec<u8>,
        codec: Codecs,
    ) {
        let arc_reader = self.inner.clone();
        let arc_deque = self.read_bufs_queue.clone();
        let send_res = self.sender.clone();
        let cur_task_num = self.cur_task;
        self.cur_task += 1;
        let cvar = self.cvar_task_num.clone();

        self.pool.spawn(move || {
            let mut read_buf = arc_deque.lock().unwrap().pop_front().unwrap();
            read_buf.resize(block_size, 0);
            let mut reader = arc_reader.lock().unwrap();
            (*reader).seek(SeekFrom::Start(seek_pos)).unwrap();
            (*reader).read_exact(&mut read_buf).unwrap();
            drop(reader); // Release mutex guard since reader is not needed anymore

            decompress_block(&read_buf[..], &mut buf[..], &codec).expect("Decompression failed.");
            arc_deque.lock().unwrap().push_back(read_buf);

            let (lock, cvar) = &*cvar;

            let mut waiting_for_task_num = lock.lock().unwrap();
            while *waiting_for_task_num != cur_task_num {
                waiting_for_task_num = cvar.wait(waiting_for_task_num).unwrap();
            }
            *waiting_for_task_num += 1;
            cvar.notify_one();

            send_res.send(buf).unwrap();
        });
    }

    /// Receives prefetched block. This is a blocking function. In case there is
    /// no uncompressed blocks in queue, the thread which called it will be
    /// blocked until uncompressed buffer comes.
    pub fn get_block(&mut self) -> Vec<u8> {
        self.receiver.recv().unwrap()
    }
}
