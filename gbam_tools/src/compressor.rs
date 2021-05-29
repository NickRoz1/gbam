use crate::SIZE_LIMIT;
use flume::{Receiver, Sender};
use rayon::ThreadPool;

use super::{Codecs, Fields};
use flate2::write::GzEncoder;
use flate2::Compression;
use lz4::EncoderBuilder;
use std::io::Write;
// use std::sync::{Arc, RwLock};

/// Accompanies compressed buffer to generate meta when written out
pub(crate) struct CompressTask {
    pub ordering_key: usize,
    pub field: Fields,
    pub num_items: u32,
    pub uncompr_size: usize,
    pub buf: Vec<u8>,
}
pub(crate) struct Compressor {
    compr_pool: ThreadPool,
    compr_data_tx: Sender<CompressTask>,
    compr_data_rx: Receiver<CompressTask>,
    /// Buffers shared among threads
    buf_tx: Sender<Vec<u8>>,
    buf_rx: Receiver<Vec<u8>>,
    // Total number of decompression queryies
    sent: usize,
    // Processed blocks number
    received: usize,
}

impl Compressor {
    pub fn new(thread_num: usize) -> Self {
        let (compr_data_tx, compr_data_rx) = flume::unbounded();
        let (buf_tx, buf_rx) = flume::unbounded();
        for _ in 0..thread_num {
            buf_tx.send(vec![0; SIZE_LIMIT]).unwrap();
            compr_data_tx
                .send(CompressTask {
                    ordering_key: 0,
                    field: Fields::RefID,
                    num_items: 0,
                    uncompr_size: 0,
                    buf: vec![0; SIZE_LIMIT],
                })
                .unwrap();
        }
        Compressor {
            compr_pool: rayon::ThreadPoolBuilder::new()
                .num_threads(thread_num)
                .build()
                .unwrap(),
            compr_data_tx,
            compr_data_rx,
            buf_tx,
            buf_rx,
            sent: 0,
            received: 0,
        }
    }

    pub fn compress_block(
        &mut self,
        ordering_key: usize,
        field: Fields,
        num_items: u32,
        uncompr_size: usize,
        data: Vec<u8>,
        codec: Codecs,
    ) {
        let buf_queue_tx = self.buf_tx.clone();
        let buf_queue_rx = self.buf_rx.clone();
        let compressed_tx = self.compr_data_tx.clone();
        self.sent += 1;
        self.compr_pool.install(|| {
            rayon::spawn(move || {
                let mut buf = buf_queue_rx.recv().unwrap();
                buf.clear();
                let compr_data = compress(&data[..uncompr_size], buf, codec);
                buf_queue_tx.send(data).unwrap();

                compressed_tx
                    .send(CompressTask {
                        ordering_key,
                        field,
                        num_items,
                        uncompr_size,
                        buf: compr_data,
                    })
                    .unwrap();
            });
        });
    }

    /// Drain completed tasks
    pub fn get_compr_block(&mut self) -> CompressTask {
        let task = self.compr_data_rx.recv().unwrap();
        // Correct for first dummy blocks
        if task.uncompr_size != 0 {
            self.received += 1;
        }
        task
    }

    /// Wait for all threads to finish and return leftovers
    pub fn finish(&mut self) -> Vec<CompressTask> {
        let mut leftovers = Vec::new();
        while self.received != self.sent {
            leftovers.push(self.get_compr_block());
        }
        leftovers
    }
}

pub fn compress(mut source: &[u8], dest: Vec<u8>, codec: Codecs) -> Vec<u8> {
    let compressed_bytes = match codec {
        Codecs::Gzip => {
            let mut encoder = GzEncoder::new(dest, Compression::default());
            encoder.write_all(source).unwrap();
            encoder.finish()
        }
        Codecs::Lz4 => {
            let default_compression: u32 = 4;
            let mut encoder = EncoderBuilder::new()
                .level(default_compression)
                .build(dest)
                .unwrap();
            std::io::copy(&mut source, &mut encoder).unwrap();
            let (_output, result) = encoder.finish();
            match result {
                Ok(()) => Ok(_output),
                Err(error) => Err(error),
            }
        }
    };
    compressed_bytes.unwrap()
}
