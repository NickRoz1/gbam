use crate::block::Block;
use crate::util::{fetch_block, inflate_data};

// This module preparses BAM blocks to parallelize decompression
use flume::{Receiver, Sender};
use rayon::spawn;
use std::cmp::{Ord, Ordering, PartialEq, PartialOrd};
use std::collections::BinaryHeap;
use std::io::Read;

#[allow(clippy::upper_case_acronyms)]
enum Status {
    Success(Block),
    EOF,
}

struct Task(usize, Status);

impl Ord for Task {
    fn cmp(&self, other: &Self) -> Ordering {
        // Smallest go first.
        other.0.cmp(&self.0)
    }
}

impl PartialOrd for Task {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        // Smallest go first.
        Some(self.cmp(other))
    }
}

impl Eq for Task {}

impl PartialEq for Task {
    fn eq(&self, other: &Self) -> bool {
        // There shouldn't be two WorkUnits with the same number in the Heap
        assert_ne!(self.0, other.0);
        false
    }
}

/// Prefetches and decompresses GBAM blocks
pub(crate) struct Readahead {
    // Decompressing threadpool.
    used_block_sender: Sender<Block>,
    ready_to_processing_rx: Receiver<Status>,
}

impl Readahead {
    pub fn new(mut thread_num: usize, mut reader: Box<dyn Read + Send + 'static>) -> Self {
        // No less than 3 threads to avoid deadlock.
        thread_num = std::cmp::max(thread_num, 3);
        let (read_bufs_send, read_bufs_recv) = flume::unbounded();
        let (used_block_sender, used_block_receiver) = flume::unbounded();
        let (completed_task_tx, sorting_blocks_rx) = flume::unbounded();
        let (ready_tasks_tx, ready_to_processing_rx) = flume::unbounded();
        for _ in 0..thread_num {
            read_bufs_send.send(Vec::new()).unwrap();
            used_block_sender.send(Block::default()).unwrap();
        }
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(thread_num)
            .build()
            .unwrap();

        // Ordering thread.
        pool.spawn(move || {
            // The heap is needed for cases when the blocks are not inflated in
            // proper order (as coming from input stream).
            let mut block_heap = BinaryHeap::<Task>::new();
            // Number of current block (ordered as read from input stream).
            let mut cur_block_num = 0;
            while let Ok(work_unit) = sorting_blocks_rx.recv() {
                block_heap.push(work_unit);
                // Fill queue with parsed blocks.
                while !block_heap.is_empty()
                    && (block_heap.peek().unwrap().0 == cur_block_num)
                    && !ready_tasks_tx.is_disconnected()
                {
                    ready_tasks_tx.send((block_heap.pop().unwrap()).1).unwrap();
                    // The block is extracted. Wait for next one.
                    cur_block_num += 1;
                }
            }
        });
        // Reading thread.
        pool.spawn(move || {
            let mut cur_task: usize = 0;
            while let Ok(mut block) = used_block_receiver.recv() {
                let mut read_buf = read_bufs_recv.recv().unwrap();
                let bytes_count = fetch_block(&mut reader, &mut read_buf, &mut block).unwrap();
                block.compressed_size = bytes_count as u64;

                let task_ready_to_sort_tx = completed_task_tx.clone();
                if bytes_count == 0 {
                    task_ready_to_sort_tx
                        .send(Task(cur_task, Status::EOF))
                        .unwrap();
                    // Reached EOF
                    return;
                }

                let read_buf_sender = read_bufs_send.clone();

                spawn(move || {
                    decompress_block(&read_buf, &mut block);
                    task_ready_to_sort_tx
                        .send(Task(cur_task, Status::Success(block)))
                        .unwrap();
                    if !read_buf_sender.is_disconnected() {
                        read_buf_sender.send(read_buf).unwrap();
                    }
                });

                cur_task += 1;
            }
        });
        Self {
            used_block_sender,
            ready_to_processing_rx,
        }
    }

    /// Receives prefetched block. This is a blocking function. In case there is
    /// no uncompressed blocks in queue, the thread which called it will be
    /// blocked until uncompressed buffer appears.
    pub fn get_block(&mut self, old_buf: Block) -> Option<Block> {
        // eprintln!("3.6.");
        if !self.used_block_sender.is_disconnected() {
            // Ignore even if it errs. Even though the check has been passed at
            // this point the threads might have been already terminated, so it
            // will err on send attempt (no available receivers).
            let _ = self.used_block_sender.send(old_buf);
        }
        // eprintln!("3.7.");
        match self.ready_to_processing_rx.recv().unwrap() {
            Status::Success(block) => Some(block),
            Status::EOF => None,
        }
        // eprintln!("3.8.");
    }
}

fn decompress_block(read_buf: &[u8], block: &mut Block) {
    let udata = block.data_mut();
    let udata_buf = udata.get_mut();
    inflate_data(read_buf, udata_buf).expect("Decompression failed.");
    udata.set_position(0);
}
