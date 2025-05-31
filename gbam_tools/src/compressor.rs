use bam_tools::record::fields::Fields;
use crate::SIZE_LIMIT;
use flume::{Receiver, Sender};
use rayon::ThreadPool;

use super::Codecs;
use flate2::write::GzEncoder;
use flate2::Compression;
// use lz4::EncoderBuilder;
use std::io::Write;

// use lz4_flex::block::{compress_into, get_maximum_output_size};
use lzzzz::lz4;

use crate::writer::BlockInfo;

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
                // Before compressing, we may encode RawSequence using 2-bit encoding.
                // Ref: https://www.biorxiv.org/content/10.1101/2025.04.08.647863v1
                let slice = &data[..block_info.uncompr_size];
                let encoded_bytes: &[u8] = if block_info.field == Fields::RawSequence {
                    match encode_sequence_2bit(slice, true, b'A') {
                        Some(words) => {
                            // Flatten u64 vector into &[u8]
                            buf.extend(words.iter().flat_map(|w| w.to_le_bytes()));
                            &buf
                        }
                        None => {
                            eprintln!("Invalid base found in RawSequence block.");
                            return;
                        }
                    }
                } else {
                    slice
                };

                let compr_data = compress(&encoded_bytes.to_vec(), buf, codec);
                // let compr_data = compress(&data[..block_info.uncompr_size], buf, codec);
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
            let mut encoder = GzEncoder::new(dest, Compression::default());
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
                Err(_) => Err(std::io::Error::new(
                    std::io::ErrorKind::Other,
                    "Compression error",
                )),
            }
        },
        Codecs::NoCompression => {
            dest.clear();
            dest.extend_from_slice(source);
            Ok(dest)
        }
    };
    compressed_bytes.unwrap()
}

/// Encode a DNA sequence (ACGT) into a vector of u64 words (2 bits per base).
/// Each u64 word contains up to 32 nucleotides.
/// Returns `None` if invalid characters are found and `allow_invalid` is false.
fn encode_sequence_2bit(seq: &[u8], allow_invalid: bool, replacement: u8) -> Option<Vec<u64>> {
    let mut result = Vec::new();
    let mut word: u64 = 0;
    let mut bits_filled = 0;

    // 2-bit mapping: A=00, C=01, G=10, T=11
    fn base_to_bits(base: u8) -> Option<u64> {
        match base {
            b'A' => Some(0b00),
            b'C' => Some(0b01),
            b'G' => Some(0b10),
            b'T' => Some(0b11),
            _ => None,
        }
    }

    for &base in seq {
        let valid_base = base_to_bits(base).or_else(|| {
            if allow_invalid {
                base_to_bits(replacement)
            } else {
                None
            }
        })?;

        word = (word << 2) | valid_base;
        bits_filled += 2;

        if bits_filled == 64 {
            result.push(word);
            word = 0;
            bits_filled = 0;
        }
    }

    if bits_filled > 0 {
        word <<= 64 - bits_filled;
        result.push(word);
    }

    Some(result)
}
