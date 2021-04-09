use crate::{u32_size, u64_size, u8_size, Fields, RawRecord, FIELDS_NUM};
use byteorder::{LittleEndian, WriteBytesExt};
use serde::{Deserialize, Serialize};
use serde_json::Result;
use std::io::{Seek, SeekFrom, Write};

static GBAM_MAGIC: &[u8] = b"geeBAM10";

#[derive(Serialize, Deserialize)]
enum CODECS {
    gzip,
    lz4,
}
#[derive(Serialize, Deserialize)]
struct BlockMeta {
    seekpos: u64,
    numitems: u32,
}

#[derive(Serialize, Deserialize)]
struct POS {
    item_size: u32,
    block_size: u32,
    codecs: CODECS,
    blocks: Vec<BlockMeta>,
}

#[derive(Serialize, Deserialize)]
struct MAPQ {
    item_size: u32,
    block_size: u32,
    codecs: CODECS,
    blocks: Vec<BlockMeta>,
}

#[derive(Serialize, Deserialize)]
struct FileMeta {
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

    fn get_blocks(&mut self, field: &Fields) -> &mut Vec<BlockMeta> {
        match field {
            Fields::Mapq => &mut self.mapq.blocks,
            Fields::Pos => &mut self.pos.blocks,
            _ => panic!("Unreachable!"),
        }
    }

    fn get_field_size(&self, field: &Fields) -> u32 {
        match field {
            Fields::Mapq => self.mapq.item_size,
            Fields::Pos => self.pos.item_size,
            _ => panic!("Unreachable!"),
        }
    }
}

const SIZE_LIMIT: usize = 16777216;
/// Groups records before writing out to file.
pub struct Writer<W>
where
    W: Write + Seek,
{
    chunks: Vec<Vec<u8>>,
    // Current item index
    offsets: [usize; FIELDS_NUM],
    num_items: [u32; FIELDS_NUM],
    fields_to_flush: [bool; FIELDS_NUM],
    file_meta: FileMeta,
    inner: W,
}

impl<W> Writer<W>
where
    W: Write + Seek,
{
    pub fn new(mut inner: W) -> Self {
        inner
            .seek(SeekFrom::Start((u64_size + u32_size * 2 + u64_size) as u64))
            .unwrap();
        Writer {
            chunks: vec![vec![0; SIZE_LIMIT]; FIELDS_NUM],
            offsets: [0; FIELDS_NUM],
            num_items: [0; FIELDS_NUM],
            fields_to_flush: [false; FIELDS_NUM],
            file_meta: FileMeta::new(),
            inner: inner,
        }
    }
    pub fn push_record(&mut self, record: &RawRecord) {
        for field in Fields::iterator().map(|v| *v as usize) {
            let cur_chunk = &mut self.chunks[field];
            let offset = &mut self.offsets[field];
            let item_counter = &mut self.num_items[field];
            let new_data = record.get_bytes(&Fields::Pos);
            cur_chunk[*offset..*offset + new_data.len()].clone_from_slice(new_data);
            *offset += new_data.len();
            *item_counter += 1;
        }

        let mut flush_required = false;
        for (idx, offset) in self.offsets.iter_mut().enumerate() {
            if *offset > SIZE_LIMIT {
                self.fields_to_flush[idx] = true;
                flush_required = true;
            }
        }

        if flush_required {
            self.flush()
        }
    }

    fn flush(&mut self) {
        let fields_to_flush = self.fields_to_flush.clone();
        for (_, field) in fields_to_flush
            .iter()
            .zip(Fields::iterator())
            .filter(|(v, field)| **v)
        {
            let meta = self.generate_meta(field);
            let field_meta = self.file_meta.get_blocks(field);
            field_meta.push(meta);
            println!("---------- {}", self.offsets[*field as usize]);
            // Write the data
            self.inner
                .write(&self.chunks[*field as usize][0..self.offsets[*field as usize]]);
        }
        for (idx, field) in self.fields_to_flush.iter_mut().enumerate() {
            if *field {
                self.offsets[idx] = 0;
                self.num_items[idx] = 0;
            }
            *field = false;
        }
    }

    fn generate_meta(&mut self, field: &Fields) -> BlockMeta {
        let item_size = self.file_meta.get_field_size(field);
        let field_meta = self.file_meta.get_blocks(field);
        let mut seek_pos = 0;
        if !field_meta.is_empty() {
            let previous = field_meta.last().unwrap();
            seek_pos = previous.seekpos + (previous.numitems * item_size) as u64;
        }
        // Space taken by values
        let offset = self.offsets[*field as usize];
        BlockMeta {
            seekpos: seek_pos,
            numitems: self.num_items[*field as usize],
        }
    }

    /// Terminates the writer. Always call after writting all the data.
    pub fn finish(&mut self) {
        // for field in self.fields_to_flush.iter_mut() {
        //     *field = true;
        // }
        self.fields_to_flush[Fields::Mapq as usize] = true;
        self.fields_to_flush[Fields::Pos as usize] = true;
        self.flush();
        let cur_pos = self.inner.seek(SeekFrom::Current(0)).unwrap();
        // Write meta
        let main_meta = serde_json::to_string(&self.file_meta).unwrap();
        self.inner.write(&main_meta.as_bytes()[..]);
        self.inner.seek(SeekFrom::Start(0)).unwrap();
        self.inner.write(GBAM_MAGIC);
        self.inner.write_u32::<LittleEndian>(1).unwrap();
        self.inner.write_u32::<LittleEndian>(0).unwrap();
        self.inner.write_u64::<LittleEndian>(cur_pos);
    }
}
