use super::meta::{BlockMeta, FileInfo, FileMeta, FILE_INFO_SIZE};
use super::SIZE_LIMIT;
use crate::{field_type, var_size_field_to_index, FieldType};
use crate::{is_data_field, u32_size, Fields, DATA_FIELDS_NUM, FIELDS_NUM};
use bit_vec::BitVec;
use byteorder::{LittleEndian, ReadBytesExt, WriteBytesExt};
use pyo3::prelude::*;
use std::io::{Read, Seek, SeekFrom};

// Represents a BAM record in which some fields may be omitted.
#[pyclass]
#[derive(Clone, Debug)]
pub struct GbamRecord {
    #[pyo3(get)]
    pub refid: Option<u32>,
    #[pyo3(get)]
    pub pos: Option<u32>,
    #[pyo3(get)]
    pub mapq: Option<u8>,
    #[pyo3(get)]
    pub bin: Option<u16>,
    #[pyo3(get)]
    pub flag: Option<u16>,
    #[pyo3(get)]
    pub next_refID: Option<u32>,
    #[pyo3(get)]
    pub next_pos: Option<u32>,
    #[pyo3(get)]
    pub read_name: Option<Vec<u8>>,
    #[pyo3(get)]
    pub cigar: Option<Vec<u32>>,
    #[pyo3(get)]
    pub seq: Option<Vec<u8>>,
    #[pyo3(get)]
    pub qual: Option<Vec<u8>>,
    #[pyo3(get)]
    pub tags: Option<Vec<u8>>,
}

impl GbamRecord {
    pub fn parse_from_bytes(&mut self, field: &Fields, mut bytes: &[u8]) {
        match field {
            Fields::RefID => self.refid = Some(bytes.read_u32::<LittleEndian>().unwrap()),
            Fields::Pos => self.pos = Some(bytes.read_u32::<LittleEndian>().unwrap()),
            Fields::Mapq => self.mapq = Some(bytes[0].to_owned()),
            Fields::Bin => self.bin = Some(bytes.read_u16::<LittleEndian>().unwrap()),
            Fields::Flags => self.flag = Some(bytes.read_u16::<LittleEndian>().unwrap()),
            Fields::NextRefID => self.next_refID = Some(bytes.read_u32::<LittleEndian>().unwrap()),
            Fields::NextPos => self.next_pos = Some(bytes.read_u32::<LittleEndian>().unwrap()),
            Fields::ReadName => self.read_name = Some(bytes.to_vec()),
            Fields::RawCigar => {
                self.cigar = Some(
                    bytes
                        .chunks(u32_size)
                        .map(|mut slice| slice.read_u32::<LittleEndian>().unwrap())
                        .collect(),
                )
            }
            Fields::RawSequence => self.seq = Some(bytes.to_vec()),
            Fields::RawQual => self.qual = Some(bytes.to_vec()),
            Fields::RawTags => self.tags = Some(bytes.to_vec()),
            _ => panic!("Not yet covered type: {}", field.to_string()),
        }
    }
}

impl Default for GbamRecord {
    fn default() -> Self {
        GbamRecord {
            refid: None,
            pos: None,
            mapq: None,
            bin: None,
            flag: None,
            next_refID: None,
            next_pos: None,
            read_name: None,
            cigar: None,
            seq: None,
            qual: None,
            tags: None,
        }
    }
}

pub trait ReadSeekSend: Read + Seek + Send {}

impl<T> ReadSeekSend for T where T: Read + Seek + Send {}
#[pyclass]
/// Reader for GBAM file format
pub struct Reader {
    inner: Box<dyn ReadSeekSend>, // Box is necessary to use with PyO3 -  no generic parameters are allowed!
    file_meta: FileMeta,
    file_info: FileInfo,
    parsing_template: ParsingTemplate,
    buffers: Vec<Vec<u8>>,
    /// Current active record
    pub cur_record: GbamRecord,
    cur_num: u64,
    cur_block: Vec<i32>,
    /// How many record have been loaded in memory already for each block. It's
    /// needed since the blocks sizes may be different.
    loaded_records_num: Vec<u64>,
}

#[pyclass]
#[derive(Clone, Debug)]
/// This struct regulates what fields are getting parsed from GBAM file.
pub struct ParsingTemplate {
    inner: Vec<Option<Fields>>,
}

impl ParsingTemplate {
    /// Create new parsing templates with all fields set to false
    pub fn new() -> Self {
        Self {
            inner: ((0..FIELDS_NUM).map(|_| None).collect()),
        }
    }
    /// Set field value
    pub fn set(&mut self, field: &Fields, val: bool) {
        match field_type(field) {
            FieldType::FixedSized => {
                self.inner[*field as usize] = Self::bool_to_val(field, val);
            }
            FieldType::VariableSized => {
                self.inner[*field as usize] = Self::bool_to_val(field, val);
                self.inner[var_size_field_to_index(field) as usize] =
                    Self::bool_to_val(&var_size_field_to_index(field), val);
            }
        }
    }

    fn bool_to_val(field: &Fields, val: bool) -> Option<Fields> {
        match val {
            true => Some(*field),
            false => None,
        }
    }
    /// Get iterator over fields currently requested for parsing
    pub fn get_active_fields_iter<'a>(&'a self) -> impl Iterator<Item = &'a Fields> {
        self.inner
            .iter()
            .filter(|x| x.is_some())
            .map(|x| x.as_ref().unwrap())
    }

    /// Get iterator over data fields (no index fields) currently requested for parsing
    pub fn get_active_data_fields_iter<'a>(&'a self) -> impl Iterator<Item = &'a Fields> {
        self.inner
            .iter()
            .filter(|x| x.is_some() && is_data_field(&x.unwrap()))
            .map(|x| x.as_ref().unwrap())
    }

    /// Get fields currently requested for parsing
    pub fn get_active_fields(&self) -> Vec<Fields> {
        self.get_active_fields_iter()
            .map(|field| field.to_owned())
            .collect::<Vec<_>>()
    }
    /// Set all fields to active state
    pub fn set_all(&mut self) {
        for (field, val) in Fields::iterator().zip(self.inner.iter_mut()) {
            if val.is_none() {
                *val = Some(*field);
            }
        }
    }
    /// Set all fields to disabled state
    pub fn clear(&mut self) {
        self.inner
            .iter_mut()
            .filter(|x| x.is_some())
            .for_each(|e| *e = None);
    }
}

#[pymethods]
impl ParsingTemplate {
    #[new]
    /// Create new ParsingTemplate with predefined values
    pub fn new_from_python(fields: Vec<bool>) -> Self {
        assert_eq!(fields.len(), DATA_FIELDS_NUM);
        let mut tmplt = ParsingTemplate::new();
        println!("{:?}", tmplt);
        for (field, val) in Fields::iterator().zip(fields) {
            println!("{:?}", field);
            tmplt.set(field, val);
        }
        tmplt
    }
}
// Methods to be called from Rust.
impl Reader {
    /// Create new reader
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
        let mut reader = Self {
            inner: inner,
            file_meta: file_meta,
            file_info: file_info,
            parsing_template: parsing_tmplt,
            buffers: (0..FIELDS_NUM).map(|_| vec![0; SIZE_LIMIT]).collect(),
            cur_record: GbamRecord::default(),
            cur_num: 0,
            cur_block: vec![-1; FIELDS_NUM], // to preload first blocks
            loaded_records_num: vec![0; FIELDS_NUM],
        };
        println!("{:?}", reader.parsing_template.get_active_fields());
        for field in reader.parsing_template.get_active_fields().iter() {
            reader.load_next_block(field);
        }
        Ok(reader)
    }

    fn load_next_block(&mut self, field: &Fields) {
        let next_block_num = self.cur_block[*field as usize] + 1;
        println!(
            "Field: {}. VALUE :{:?}",
            field.to_string(),
            self.cur_block[*field as usize]
        );
        let block_size = match field_type(field) {
            FieldType::VariableSized => {
                self.file_meta.get_blocks_sizes(field)[next_block_num as usize]
            }
            FieldType::FixedSized => {
                let item_size = self.file_meta.get_field_size(field).unwrap();
                let next_block_meta = &self.file_meta.get_blocks(field)[next_block_num as usize];
                next_block_meta.numitems * item_size
            }
        };

        let seekpos = self.file_meta.get_blocks(field)[next_block_num as usize].seekpos;
        self.load_data(seekpos, field, block_size).unwrap();

        self.cur_block[*field as usize] = next_block_num;
        println!("{:?}", self.cur_block[*field as usize]);
        self.loaded_records_num[*field as usize] +=
            self.file_meta.get_blocks(field)[next_block_num as usize].numitems as u64;
    }

    fn load_data(&mut self, seekpos: u64, field: &Fields, block_size: u32) -> std::io::Result<()> {
        self.inner.seek(SeekFrom::Start(seekpos))?;
        let buf = &mut self.buffers[*field as usize];
        buf.resize(block_size as usize, 0);
        self.inner.read_exact(buf)?;
        Ok(())
    }

    /// Calculates offset to the record inside the block.
    fn get_offset_into_block(&self, field: &Fields) -> u32 {
        let rec_num_in_block = self.calc_rec_num_in_block(field);
        match field_type(field) {
            FieldType::VariableSized => {
                if rec_num_in_block == 0 {
                    0
                } else {
                    let idx_field = &var_size_field_to_index(field);
                    let idx_field_rec_num = self.calc_rec_num_in_block(idx_field);
                    self.get_idx_val(idx_field, idx_field_rec_num - 1) as u32
                }
            }
            FieldType::FixedSized => {
                rec_num_in_block as u32 * self.file_meta.get_field_size(field).unwrap() as u32
            }
        }
    }

    fn calc_rec_num_in_block(&self, field: &Fields) -> u64 {
        println!("{:?}", field);
        assert!(self.cur_block[*field as usize] > -1);
        let block_meta =
            &self.file_meta.view_blocks(field)[self.cur_block[*field as usize] as usize];
        // self.cur_num is universal, but amount of items in buffer is not
        let dist_from_back_of_block = self.loaded_records_num[*field as usize] - self.cur_num;
        let rec_num_in_block = block_meta.numitems as u64 - dist_from_back_of_block;
        rec_num_in_block
    }

    fn is_last_in_block(&self, rec_num: u64, field: &Fields) -> bool {
        assert!(self.cur_block[*field as usize] > -1);
        let block_meta =
            &self.file_meta.view_blocks(field)[self.cur_block[*field as usize] as usize];
        (block_meta.numitems - 1) as u64 == rec_num
    }

    /// Parses current index value. All these fields are u32.
    fn get_idx_val(&self, field: &Fields, rec_num_in_block: u64) -> u32 {
        println!("{:?}", field);
        match field {
            Fields::LName
            | Fields::SequenceLength
            | Fields::RawSeqLen
            | Fields::RawTagsLen
            | Fields::NCigar => self
                .get_field_bytes(field, rec_num_in_block)
                .read_u32::<LittleEndian>()
                .unwrap(),
            _ => panic!("Only index fields allowed!"),
        }
    }

    fn get_field_bytes(&self, field: &Fields, rec_num_in_block: u64) -> &[u8] {
        let offset = self.get_offset_into_block(field) as usize;
        let buffer = &self.buffers[*field as usize];
        let item_size = self.get_cur_item_size(field, rec_num_in_block) as usize;
        &buffer[offset..offset + item_size]
    }

    /// Size of current item in bytes
    fn get_cur_item_size(&self, field: &Fields, rec_num_in_block: u64) -> usize {
        match field_type(field) {
            FieldType::VariableSized => {
                let idx_field = var_size_field_to_index(field);
                let cur_idx_rec_num = self.calc_rec_num_in_block(&idx_field);
                if rec_num_in_block == 0 {
                    self.get_idx_val(&idx_field, cur_idx_rec_num) as usize
                } else {
                    let item_size = self.get_idx_val(&idx_field, cur_idx_rec_num)
                        - self.get_idx_val(&idx_field, cur_idx_rec_num - 1);
                    item_size as usize
                }
            }
            FieldType::FixedSized => self.file_meta.get_field_size(field).unwrap() as usize,
        }
    }
    fn block_exhausted(&self, field: &Fields) -> bool {
        self.cur_num == self.loaded_records_num[*field as usize]
    }

    /// Get next GBAM record
    pub fn next(&mut self) -> Option<&GbamRecord> {
        for field in self.parsing_template.get_active_fields() {
            let cur_block_num = self.cur_block[field as usize];
            if self.block_exhausted(&field) {
                if (cur_block_num + 1) as usize == self.file_meta.get_blocks(&field).len() {
                    // No more records of this type in buffer
                    return None;
                }
                // One of the blocks have been exhausted;
                self.load_next_block(&field);
            }
        }

        for active_field in self.parsing_template.get_active_data_fields_iter() {
            let offset = self.get_offset_into_block(active_field) as usize;
            let buffer = &self.buffers[*active_field as usize];
            let rec_num_in_block = self.calc_rec_num_in_block(active_field);
            let item_size = self.get_cur_item_size(active_field, rec_num_in_block) as usize;

            println!("FIELD {} : {:?}", active_field.to_string(), buffer.len());
            self.cur_record
                .parse_from_bytes(active_field, &buffer[offset..offset + item_size]);
        }
        // TODO :: NULLIFY CURNUM ON REPEATED CALLS
        self.cur_num += 1;
        Some(&self.cur_record)
    }
}

// // Methods to be called from Python.
#[pymethods]
impl Reader {
    fn next_rec(&mut self) -> PyResult<Option<GbamRecord>> {
        let next = self.next().cloned();
        Ok(next)
    }

    #[new]
    /// Create new reader for file at path
    pub fn new_for_file(path: &str, tmplt: ParsingTemplate) -> Self {
        use std::fs::File;
        use std::io::{self, prelude::*, BufReader};

        let file = File::open(path).unwrap();
        let reader = BufReader::new(file);
        Self::new(Box::new(reader), tmplt).unwrap()
    }
}
