use super::meta::{BlockMeta, FileInfo, FileMeta, FILE_INFO_SIZE};
use super::SIZE_LIMIT;
use crate::{Fields, FIELDS_NUM};
use bit_vec::BitVec;
use byteorder::{LittleEndian, ReadBytesExt, WriteBytesExt};
use pyo3::prelude::*;
use std::io::{Read, Seek, SeekFrom};

// Represents a BAM record in which some fields may be omitted.
#[pyclass]
pub struct GbamRecord {
    pub mapq: Option<u8>,
    pub pos: Option<u32>,
}

impl GbamRecord {
    pub fn parse_from_bytes(&mut self, field: &Fields, mut bytes: &[u8]) {
        match field {
            Fields::Mapq => self.mapq = Some(bytes[0].to_owned()),
            Fields::Pos => self.pos = Some(bytes.read_u32::<LittleEndian>().unwrap()),
            _ => panic!("Not yet covered type."),
        }
    }
}

pub trait ReadSeekSend: Read + Seek + Send {}

impl<T> ReadSeekSend for T where T: Read + Seek + Send {}
#[pyclass]
pub struct Reader {
    inner: Box<dyn ReadSeekSend>, // Necessary to use with PyO3. No generic parameters are allowed!
    file_meta: FileMeta,
    file_info: FileInfo,
    parsing_template: ParsingTemplate,
    buffers: Vec<Vec<u8>>,
    cur_record: GbamRecord,
    cur_num: u64,
    cur_block: Vec<u64>,
    /// How many record have been loaded in memory already for each block. It's
    /// needed since the blocks sizes may be different.
    loaded_records_num: Vec<u64>,
}

// This vector regulates what fields are getting parsed from GBAM file.
pub struct ParsingTemplate(Vec<Option<Fields>>);

impl ParsingTemplate {
    pub fn new() -> Self {
        Self((0..FIELDS_NUM).map(|_| None).collect())
    }
    pub fn set(&mut self, field: &Fields, val: bool) {
        self.0[*field as usize] = match val {
            true => Some(*field),
            false => None,
        }
    }

    pub fn get_active_fields<'a>(&'a self) -> impl Iterator<Item = &'a Fields> {
        self.0
            .iter()
            .filter(|x| x.is_some())
            .map(|x| x.as_ref().unwrap())
    }

    pub fn set_all(&mut self) {
        for (field, val) in Fields::iterator().zip(self.0.iter_mut()) {
            if val.is_none() {
                *val = Some(*field);
            }
        }
    }

    pub fn clear(&mut self) {
        self.0
            .iter_mut()
            .filter(|x| x.is_some())
            .for_each(|e| *e = None);
    }
}
// Methods to be called from Rust.
impl Reader {
    pub fn new(
        mut inner: Box<dyn ReadSeekSend>,
        parsing_tmplt: ParsingTemplate,
    ) -> std::io::Result<Self> {
        let mut file_info_bytes: [u8; FILE_INFO_SIZE] = [0; FILE_INFO_SIZE];
        inner.read_exact(&mut file_info_bytes[..])?;
        let file_info = FileInfo::from(&file_info_bytes[..]);
        // Read file meta
        inner.seek(SeekFrom::Start(file_info.seekpos))?;
        let mut buf = Vec::<u8>::new();
        inner.read_to_end(&mut buf)?;

        let file_meta_json_str = String::from_utf8(buf).unwrap();
        let file_meta =
            serde_json::from_str(&file_meta_json_str).expect("File meta json string was damaged.");
        Ok(Self {
            inner: inner,
            file_meta: file_meta,
            file_info: file_info,
            parsing_template: parsing_tmplt,
            buffers: (0..FIELDS_NUM).map(|_| vec![0; SIZE_LIMIT]).collect(),
            cur_record: GbamRecord {
                mapq: None,
                pos: None,
            },
            cur_num: 0,
            cur_block: vec![0; FIELDS_NUM],
            loaded_records_num: vec![0; FIELDS_NUM],
        })
    }

    fn load_next_block(&mut self, field: &Fields) {
        let next_block_num = self.cur_block[*field as usize] + 1;
        let item_size = self.file_meta.get_field_size(field);
        let next_block_meta = &self.file_meta.get_blocks(field)[next_block_num as usize];
        let block_size = next_block_meta.numitems * item_size;

        self.inner.seek(SeekFrom::Start(next_block_meta.seekpos));
        let buf = &mut self.buffers[*field as usize];
        buf.resize(block_size as usize, 0);
        self.inner.read_exact(buf);
    }

    /// Calculates offset to the record inside the block.
    fn get_offset_into_block(&self, field: &Fields) -> u64 {
        let block_meta =
            &self.file_meta.view_blocks(field)[self.cur_block[*field as usize] as usize];

        let dist_from_back_of_block = self.loaded_records_num[*field as usize] - self.cur_num;
        let rec_num_in_block = block_meta.numitems as u64 - dist_from_back_of_block;

        let offset = rec_num_in_block as u64 * self.file_meta.get_field_size(field) as u64;
        offset
    }

    pub fn next(&mut self) -> Option<&GbamRecord> {
        for field in Fields::iterator().filter(|v| **v == Fields::Mapq || **v == Fields::Pos) {
            let cur_block_num = self.cur_block[*field as usize];
            if self.cur_num == self.loaded_records_num[*field as usize] {
                if (cur_block_num + 1) as usize == self.file_meta.get_blocks(field).len() {
                    return None;
                }
                // One of the blocks have been exhausted;
                self.load_next_block(field);
            }
        }

        for active_field in self.parsing_template.get_active_fields() {
            let offset = self.get_offset_into_block(active_field) as usize;
            let buffer = &self.buffers[*active_field as usize];
            let item_size = self.file_meta.get_field_size(active_field) as usize;
            self.cur_record
                .parse_from_bytes(active_field, &buffer[offset..offset + item_size]);
        }

        Some(&self.cur_record)
    }
}

// Methods to be called from Python.
#[pymethods]
impl Reader {}
