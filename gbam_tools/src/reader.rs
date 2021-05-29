use super::meta::{Codecs, FileInfo, FileMeta, FILE_INFO_SIZE};
use super::writer::calc_crc_for_meta_bytes;
#[cfg(feature = "python-ffi")]
use super::DATA_FIELDS_NUM;
use super::SIZE_LIMIT;
use super::{decode_cigar, decode_seq, is_data_field, Fields, FIELDS_NUM};
use super::{field_type, var_size_field_to_index, FieldType};
use byteorder::{LittleEndian, ReadBytesExt};
use flate2::write::GzDecoder;
use lz4::Decoder;
#[cfg(feature = "python-ffi")]
use pyo3::prelude::*;
use std::io::{Read, Seek, SeekFrom};

#[cfg(not(feature = "python-ffi"))]
#[derive(Clone, Debug)]
/// Represents a GBAM record in which some fields may be omitted.
pub struct GbamRecord {
    /// Reference sequence ID
    pub refid: Option<u32>,
    /// 0-based leftmost coordinate
    pub pos: Option<u32>,
    /// Mapping quality
    pub mapq: Option<u8>,
    /// BAI index bin,
    pub bin: Option<u16>,
    /// Bitwise flags
    pub flag: Option<u16>,
    /// Ref-ID of the next segment
    pub next_ref_id: Option<i32>,
    /// 0-based leftmost pos of the next segmen
    pub next_pos: Option<i32>,
    /// Template length
    pub tlen: Option<i32>,
    /// Read name
    pub read_name: Option<Vec<u8>>,
    /// CIGAR
    pub cigar: Option<String>,
    /// 4-bit  encoded  read
    pub seq: Option<String>,
    /// Phred-scaled base qualities.
    pub qual: Option<Vec<u8>>,
    /// List of auxiliary data
    pub tags: Option<Vec<u8>>,
}

/// cfg_attribute currently doesn't work with pyo3(get)
#[cfg(feature = "python-ffi")]
#[pyclass]
#[derive(Clone, Debug)]
/// Represents a GBAM record in which some fields may be omitted.
pub struct GbamRecord {
    /// Reference sequence ID
    #[pyo3(get)]
    pub refid: Option<u32>,
    /// 0-based leftmost coordinate
    #[pyo3(get)]
    pub pos: Option<u32>,
    /// Mapping quality
    #[pyo3(get)]
    pub mapq: Option<u8>,
    /// BAI index bin,
    #[pyo3(get)]
    pub bin: Option<u16>,
    /// Bitwise flags
    #[pyo3(get)]
    pub flag: Option<u16>,
    /// Ref-ID of the next segment
    #[pyo3(get)]
    pub next_ref_id: Option<i32>,
    /// 0-based leftmost pos of the next segmen
    #[pyo3(get)]
    pub next_pos: Option<i32>,
    /// Template length
    #[pyo3(get)]
    pub tlen: Option<i32>,
    /// Read name
    #[pyo3(get)]
    pub read_name: Option<Vec<u8>>,
    /// CIGAR
    #[pyo3(get)]
    pub cigar: Option<String>,
    /// 4-bit  encoded  read
    #[pyo3(get)]
    pub seq: Option<String>,
    /// Phred-scaled base qualities.
    #[pyo3(get)]
    pub qual: Option<Vec<u8>>,
    /// List of auxiliary data
    #[pyo3(get)]
    pub tags: Option<Vec<u8>>,
}

// TODO :: ADD TEMPLATE LENGTHS TO GBAM RECORD
// TODO :: REMOVE CG TAG FROM ORIGINAL FILE
impl GbamRecord {
    pub(crate) fn parse_from_bytes(&mut self, field: &Fields, mut bytes: &[u8]) {
        match field {
            Fields::RefID => self.refid = Some(bytes.read_u32::<LittleEndian>().unwrap()),
            Fields::Pos => self.pos = Some(bytes.read_u32::<LittleEndian>().unwrap()),
            Fields::Mapq => self.mapq = Some(bytes[0].to_owned()),
            Fields::Bin => self.bin = Some(bytes.read_u16::<LittleEndian>().unwrap()),
            Fields::Flags => self.flag = Some(bytes.read_u16::<LittleEndian>().unwrap()),
            Fields::NextRefID => self.next_ref_id = Some(bytes.read_i32::<LittleEndian>().unwrap()),
            Fields::NextPos => self.next_pos = Some(bytes.read_i32::<LittleEndian>().unwrap()),
            Fields::TemplateLength => self.tlen = Some(bytes.read_i32::<LittleEndian>().unwrap()),
            Fields::ReadName => self.read_name = Some(bytes.to_vec()),
            Fields::RawCigar => {
                // self.cigar = Some(
                //     bytes
                //         .chunks(U32_SIZE)
                //         .map(|mut slice| slice.read_u32::<LittleEndian>().unwrap())
                //         .collect(),
                // )
                self.cigar = Some(decode_cigar(bytes))
            }
            Fields::RawSequence => self.seq = Some(decode_seq(bytes)),
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
            next_ref_id: None,
            next_pos: None,
            tlen: None,
            read_name: None,
            cigar: None,
            seq: None,
            qual: None,
            tags: None,
        }
    }
}

impl std::fmt::Display for GbamRecord {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        std::fmt::Debug::fmt(self, f)
    }
}

#[cfg_attr(feature = "python-ffi", pymethods)]
impl GbamRecord {
    /// Convert GBAM structure to string representation
    pub fn to_str(&self) -> String {
        self.to_string()
    }
}

/// These trait has to be implemented for any inner reader supplying data to GBAM reader
pub trait ReadSeekSend: Read + Seek + Send {}

impl<T> ReadSeekSend for T where T: Read + Seek + Send {}
/// Reader for GBAM file format
#[cfg_attr(feature = "python-ffi", pyclass)]
pub struct Reader {
    inner: Box<dyn ReadSeekSend>, // Box is necessary to use with PyO3 -  no generic parameters are allowed!
    file_meta: FileMeta,
    parsing_template: ParsingTemplate,
    buffers: Vec<Vec<u8>>,
    // Decompression buffer to avoid allocations
    decompr_buf: Vec<u8>,
    /// Current active record
    pub cur_record: GbamRecord,
    cur_num: u64,
    cur_block: Vec<i32>,
    /// How many record have been loaded in memory already for each block. It's
    /// needed since the blocks sizes may be different.
    loaded_records_num: Vec<u64>,
    loaded_before_this_block: Vec<u64>,
    /// When parsing variable sized record, two indexes are needed. However, the
    /// indexes are splitted two, creating a situation when the previous index
    /// was already flushed with it's block, but it's required required to parse
    /// a variable sized field. To avoid that, last index field of each block is
    /// saved.
    prev_block_last_idx: Option<u32>,
}

#[cfg_attr(feature = "python-ffi", pyclass)]
/// This struct regulates what fields are getting parsed from GBAM file.

#[derive(Clone, Debug)]
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
    #[allow(clippy::needless_lifetimes)]
    pub fn get_active_fields_iter<'a>(&'a self) -> impl Iterator<Item = &'a Fields> {
        self.inner
            .iter()
            .filter(|x| x.is_some())
            .map(|x| x.as_ref().unwrap())
    }

    /// Get iterator over data fields (no index fields) currently requested for parsing
    #[allow(clippy::needless_lifetimes)]
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

impl Default for ParsingTemplate {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(feature = "python-ffi")]
#[pymethods]
impl ParsingTemplate {
    #[new]
    /// Create new ParsingTemplate with predefined values
    pub fn new_from_python(fields: Vec<bool>) -> Self {
        assert_eq!(fields.len(), DATA_FIELDS_NUM);
        let mut tmplt = ParsingTemplate::new();
        for (field, val) in Fields::iterator().zip(fields) {
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
        if calc_crc_for_meta_bytes(&buf[..]) != file_info.crc32 {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                "Metadata JSON was damaged.",
            ));
        }
        let file_meta_json_str = String::from_utf8(buf).unwrap();
        let file_meta =
            serde_json::from_str(&file_meta_json_str).expect("File meta json string was damaged.");
        let mut reader = Self {
            inner,
            file_meta,
            parsing_template: parsing_tmplt,
            buffers: (0..FIELDS_NUM).map(|_| vec![0; SIZE_LIMIT]).collect(),
            decompr_buf: Vec::<u8>::new(),
            cur_record: GbamRecord::default(),
            cur_num: 0,
            cur_block: vec![-1; FIELDS_NUM], // to preload first blocks
            loaded_records_num: vec![0; FIELDS_NUM],
            loaded_before_this_block: vec![0; FIELDS_NUM],
            prev_block_last_idx: None,
        };

        for field in reader.parsing_template.get_active_fields().iter() {
            reader.load_next_block(field);
        }
        Ok(reader)
    }

    /// Create new reader for file at path
    pub fn new_for_file(path: &str, tmplt: ParsingTemplate) -> Self {
        use std::fs::File;
        use std::io::BufReader;

        let file = File::open(path).unwrap();
        let reader = BufReader::new(file);
        Self::new(Box::new(reader), tmplt).unwrap()
    }

    fn load_next_block(&mut self, field: &Fields) {
        let next_block_num = self.cur_block[*field as usize] + 1;
        let block_size = self.file_meta.get_blocks_sizes(field)[next_block_num as usize];

        let seekpos = self.file_meta.get_blocks(field)[next_block_num as usize].seekpos;
        self.load_data(seekpos, field, block_size).unwrap();

        self.cur_block[*field as usize] = next_block_num;

        self.loaded_before_this_block[*field as usize] = self.loaded_records_num[*field as usize];
        self.loaded_records_num[*field as usize] +=
            self.file_meta.get_blocks(field)[next_block_num as usize].numitems as u64;
    }

    fn load_data(&mut self, seekpos: u64, field: &Fields, block_size: u32) -> std::io::Result<()> {
        self.inner.seek(SeekFrom::Start(seekpos))?;
        self.decompr_buf.resize(block_size as usize, 0);
        self.inner.read_exact(&mut self.decompr_buf)?;
        Self::decompress_block(
            &self.decompr_buf,
            &mut self.buffers[*field as usize],
            self.file_meta.get_field_codec(field),
        )?;
        Ok(())
    }

    fn decompress_block(source: &[u8], mut dest: &mut [u8], codec: &Codecs) -> std::io::Result<()> {
        use std::io::Write;
        match codec {
            Codecs::Gzip => {
                let mut decoder = GzDecoder::new(dest);
                decoder.write_all(source)?;
                decoder.try_finish()?;
            }
            Codecs::Lz4 => {
                let mut decoder = Decoder::new(source)?;
                std::io::copy(&mut decoder, &mut dest)?;
            }
        };
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
                    if idx_field_rec_num == 0 {
                        // We are on the split of blocks
                        self.prev_block_last_idx.unwrap()
                    } else {
                        self.get_idx_val(idx_field, idx_field_rec_num - 1) as u32
                    }
                }
            }
            FieldType::FixedSized => {
                rec_num_in_block as u32 * self.file_meta.get_field_size(field).unwrap() as u32
            }
        }
    }

    fn calc_rec_num_in_block(&self, field: &Fields) -> u64 {
        assert!(self.cur_block[*field as usize] > -1);
        // let block_meta =
        // &self.file_meta.view_blocks(field)[self.cur_block[*field as usize] as usize];
        // self.cur_num is universal, but amount of items in buffer is not
        // let dist_from_back_of_block = self.loaded_records_num[*field as usize] - self.cur_num;
        self.cur_num - self.loaded_before_this_block[*field as usize]
        // block_meta.numitems as u64 - dist_from_back_of_block
    }

    #[allow(dead_code)]
    fn is_last_in_block(&self, rec_num: u64, field: &Fields) -> bool {
        assert!(self.cur_block[*field as usize] > -1);
        let block_meta =
            &self.file_meta.view_blocks(field)[self.cur_block[*field as usize] as usize];
        (block_meta.numitems - 1) as u64 == rec_num
    }

    /// Parses current index value. All these fields are u32.
    fn get_idx_val(&self, field: &Fields, rec_num_in_block: u64) -> u32 {
        if is_data_field(field) {
            panic!("Only index fields allowed!")
        };
        self.get_idx_field_bytes_for_rec(field, rec_num_in_block)
            .read_u32::<LittleEndian>()
            .unwrap()
    }

    /// Query bytes for index field containing value for specified record number
    fn get_idx_field_bytes_for_rec(&self, field: &Fields, rec_num_in_block: u64) -> &[u8] {
        if is_data_field(field) {
            panic!("Only index fields allowed!")
        }
        let item_size = self.get_cur_item_size(field, rec_num_in_block) as usize;
        let offset = item_size * rec_num_in_block as usize;
        let buffer = &self.buffers[*field as usize];
        &buffer[offset..offset + item_size]
    }

    #[allow(dead_code)]
    fn get_variable_field_bytes(&self, field: &Fields, rec_num_in_block: u64) -> &[u8] {
        let offset = self.get_offset_into_block(field) as usize;
        let buffer = &self.buffers[*field as usize];
        let item_size = self.get_cur_item_size(field, rec_num_in_block) as usize;
        &buffer[offset..offset + item_size]
    }

    /// Size of current item in bytes
    /// TODO: Separate function for constant size fields
    fn get_cur_item_size(&self, field: &Fields, rec_num_in_block: u64) -> usize {
        match field_type(field) {
            FieldType::VariableSized => {
                let idx_field = var_size_field_to_index(field);
                let cur_idx_rec_num = self.calc_rec_num_in_block(&idx_field);
                if rec_num_in_block == 0 {
                    self.get_idx_val(&idx_field, cur_idx_rec_num) as usize
                } else {
                    let cur_end = self.get_idx_val(&idx_field, cur_idx_rec_num);
                    let prev_end = match cur_idx_rec_num {
                        // We are at the block split
                        0 => self.prev_block_last_idx.unwrap(),
                        _ => self.get_idx_val(&idx_field, cur_idx_rec_num - 1),
                    };
                    let item_size = cur_end - prev_end;
                    item_size as usize
                }
            }
            FieldType::FixedSized => self.file_meta.get_field_size(field).unwrap() as usize,
        }
    }
    fn block_exhausted(&mut self, field: &Fields) -> bool {
        self.cur_num == self.loaded_records_num[*field as usize]
    }

    /// Get next GBAM record
    pub fn next_rec(&mut self) -> Option<&GbamRecord> {
        for field in self.parsing_template.get_active_fields() {
            if self.block_exhausted(&field) {
                let cur_block_num = self.cur_block[field as usize];
                if (cur_block_num + 1) as usize == self.file_meta.get_blocks(&field).len() {
                    // No more records of this type in input
                    return None;
                }
                // Save last index for next calculations.
                #[allow(clippy::bool_comparison)]
                if is_data_field(&field) == false {
                    // It's exhausted, so last rec num would be current - 1.
                    let last_idx = self.calc_rec_num_in_block(&field) - 1;
                    self.prev_block_last_idx = Some(self.get_idx_val(&field, last_idx));
                }
                self.load_next_block(&field);
            }
        }

        for active_field in self.parsing_template.get_active_data_fields_iter() {
            let offset = self.get_offset_into_block(active_field) as usize;
            let buffer = &self.buffers[*active_field as usize];
            let rec_num_in_block = self.calc_rec_num_in_block(active_field);
            let item_size = self.get_cur_item_size(active_field, rec_num_in_block) as usize;

            self.cur_record
                .parse_from_bytes(active_field, &buffer[offset..offset + item_size]);
        }
        // TODO :: NULLIFY CURNUM ON REPEATED CALLS IF WE REUSE THE READER
        self.cur_num += 1;
        Some(&self.cur_record)
    }
}

// Methods to be called from Python.
#[cfg(feature = "python-ffi")]
#[pymethods]
impl Reader {
    #[allow(clippy::unnecessary_wraps)]
    fn next_record(&mut self) -> PyResult<Option<GbamRecord>> {
        let next = self.next_rec().cloned();
        Ok(next)
    }

    /// Create new reader for file at path
    #[new]
    pub fn new_reader(path: &str, tmplt: ParsingTemplate) -> Self {
        Self::new_for_file(path, tmplt)
    }
}
