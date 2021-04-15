use super::meta::{FileInfo, FileMeta, FILE_INFO_SIZE};
use super::SIZE_LIMIT;
use crate::{Fields, FIELDS_NUM};
use bit_vec::BitVec;
use pyo3::prelude::*;
use std::io::{Read, Seek, SeekFrom};

// Represents a BAM record in which some fields may be omitted.
#[pyclass]
pub struct GbamRecord {
    mapq: Option<u8>,
    pos: Option<u32>,
}

pub trait ReadSeekSend: Read + Seek + Send {}
#[pyclass]
pub struct Reader {
    inner: Box<dyn ReadSeekSend>, // Necessary to use with PyO3. No generic parameters are allowed!
    file_meta: FileMeta,
    file_info: FileInfo,
    parsing_template: ParsingTemplate,
    buffers: Vec<Vec<u8>>,
    cur_num: u64,
    cur_record: GbamRecord,
}

// This vector regulates what fields are getting parsed from GBAM file.
pub struct ParsingTemplate(BitVec);

impl ParsingTemplate {
    pub fn new() -> Self {
        Self(BitVec::from_elem(FIELDS_NUM, false))
    }
    pub fn set(&mut self, field: &Fields, val: bool) {
        self.0.set(*field as usize, val);
    }

    pub fn clear(&mut self) {
        self.0.clear();
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
        })
    }

    fn load_next_blocks(&mut self) {}
}

// Methods to be called from Python.
#[pymethods]
impl Reader {}
