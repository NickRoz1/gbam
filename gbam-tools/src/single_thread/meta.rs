use super::SIZE_LIMIT;
use crate::GBAM_MAGIC;
use crate::{u32_size, u64_size, u8_size, Fields};
use byteorder::{LittleEndian, ReadBytesExt, WriteBytesExt};
use serde::{Deserialize, Serialize};
use serde_json::Result;
use std::io::{Read, Write};

/// Holds data related to GBAM file: gbam version, seekpos to meta.
pub struct FileInfo {
    pub gbam_version: [u32; 2],
    pub seekpos: u64,
}

impl FileInfo {
    pub fn new(gbam_version: [u32; 2], seekpos: u64) -> Self {
        FileInfo {
            gbam_version: gbam_version,
            seekpos: seekpos,
        }
    }
}

/// The GBAM magic size is 8 bytes (u64_size).
pub const FILE_INFO_SIZE: usize = u64_size + u64_size + u32_size * 2;

impl From<&[u8]> for FileInfo {
    fn from(mut bytes: &[u8]) -> Self {
        assert!(
            bytes.len() >= FILE_INFO_SIZE,
            "Not enough bytes to form file info struct.",
        );
        assert_eq!(&bytes[..u64_size], &GBAM_MAGIC[..]);
        FileInfo {
            gbam_version: [
                bytes
                    .read_u32::<LittleEndian>()
                    .expect("file info is damaged: unable to read GBAM version."),
                bytes
                    .read_u32::<LittleEndian>()
                    .expect("file info is damaged: unable to read GBAM version."),
            ],
            seekpos: bytes
                .read_u64::<LittleEndian>()
                .expect("file info is damaged: unable to read seekpos."),
        }
    }
}

impl Into<Vec<u8>> for FileInfo {
    fn into(self) -> Vec<u8> {
        let mut res = Vec::<u8>::new();
        res.write(&GBAM_MAGIC[..]);
        for val in self.gbam_version.into_iter() {
            res.write_u32::<LittleEndian>(*val);
        }
        res.write_u64::<LittleEndian>(self.seekpos).unwrap();
        res
    }
}

#[derive(Serialize, Deserialize)]
pub enum CODECS {
    gzip,
    lz4,
}
#[derive(Serialize, Deserialize)]
pub struct BlockMeta {
    pub seekpos: u64,
    pub numitems: u32,
}

#[derive(Serialize, Deserialize)]
pub struct POS {
    item_size: u32,
    block_size: u32,
    codecs: CODECS,
    blocks: Vec<BlockMeta>,
}

#[derive(Serialize, Deserialize)]
pub struct MAPQ {
    item_size: u32,
    block_size: u32,
    codecs: CODECS,
    blocks: Vec<BlockMeta>,
}

#[derive(Serialize, Deserialize)]
pub struct FileMeta {
    pos: POS,
    mapq: MAPQ,
}

impl FileMeta {
    pub fn new() -> Self {
        FileMeta {
            pos: POS {
                item_size: u32_size as u32,
                block_size: SIZE_LIMIT as u32,
                codecs: CODECS::gzip,
                blocks: Vec::<BlockMeta>::new(),
            },
            mapq: MAPQ {
                item_size: u8_size as u32,
                block_size: SIZE_LIMIT as u32,
                codecs: CODECS::gzip,
                blocks: Vec::<BlockMeta>::new(),
            },
        }
    }

    /// Used to retrieve BlockMeta vector mutable borrow, to push new blocks
    /// directly into it, avoiding field matching.
    pub fn get_blocks(&mut self, field: &Fields) -> &mut Vec<BlockMeta> {
        match field {
            Fields::Mapq => &mut self.mapq.blocks,
            Fields::Pos => &mut self.pos.blocks,
            _ => panic!("Unreachable!"),
        }
    }

    pub fn view_blocks(&self, field: &Fields) -> &Vec<BlockMeta> {
        match field {
            Fields::Mapq => &self.mapq.blocks,
            Fields::Pos => &self.pos.blocks,
            _ => panic!("Unreachable!"),
        }
    }

    pub fn get_field_size(&self, field: &Fields) -> u32 {
        match field {
            Fields::Mapq => self.mapq.item_size,
            Fields::Pos => self.pos.item_size,
            _ => panic!("Unreachable!"),
        }
    }
}
