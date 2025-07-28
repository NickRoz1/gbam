use super::Codecs;
use crate::writer::BlockInfo;
use crate::SIZE_LIMIT;
use flume::{Receiver, Sender};
use rayon::ThreadPool;

use brotli::CompressorWriter;
use flate2::write::GzEncoder;
use flate2::Compression;
use lzzzz::lz4;
use xz2::write::XzEncoder;
use zstd::stream::encode_all;

use std::io::Write;

pub(crate) enum OrderingKey {
    Key(u64),
    UnusedBlock,
}

/// Accompanies compressed buffer to generate meta when written out
pub(crate) struct CompressTask {
    pub ordering_key: OrderingKey,
    pub block_info: BlockInfo,
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
                    ordering_key: OrderingKey::UnusedBlock,
                    block_info: BlockInfo::default(),
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
        ordering_key: OrderingKey,
        block_info: BlockInfo,
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
                let compr_data = compress(&data[..block_info.uncompr_size], buf, codec);
                buf_queue_tx.send(data).unwrap();

                compressed_tx
                    .send(CompressTask {
                        ordering_key,
                        block_info,
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
        if let OrderingKey::Key(_) = task.ordering_key {
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

pub fn compress(source: &[u8], mut dest: Vec<u8>, codec: Codecs) -> Vec<u8> {
    let compressed_bytes = match codec {
        Codecs::Gzip => {
            let mut encoder = GzEncoder::new(dest, Compression::new(9));
            encoder.write_all(source).unwrap();
            encoder.finish()
        }
        Codecs::Lz4 => {
            dest.clear();
            let res = lz4::compress_to_vec(source, &mut dest, lz4::ACC_LEVEL_DEFAULT);
            match res {
                Ok(size) => {
                    dest.resize(size, 0);
                    Ok(dest)
                }
                Err(_) => Err(std::io::Error::other(
                    "Compression error",
                )),
            }
        }
        Codecs::Brotli => {
            dest.clear();
            {
                let mut writer = CompressorWriter::new(&mut dest, 4096, 6, 22);
                writer.write_all(source).unwrap();
                writer.flush().unwrap();
            }
            Ok(dest)
        }
        Codecs::Xz => {
            let mut encoder = XzEncoder::new(Vec::new(), 6);
            encoder.write_all(source).unwrap();
            let compressed = encoder.finish().unwrap();
            Ok(compressed)
        }
        Codecs::Zstd => {
            // encode_all returns a Vec<u8>
            match encode_all(source, 14) {
                Ok(c) => Ok(c),
                Err(_) => Err(std::io::Error::other(
                    "Zstd compression error",
                )),
            }
        }
        Codecs::NoCompression => {
            dest.clear();
            dest.extend_from_slice(source);
            Ok(dest)
        }
    };
    compressed_bytes.unwrap()
}
