#[cfg(feature = "python-ffi")]
use pyo3::prelude::*;

use bam_tools::record::{
    bamrawrecord::{decode_cigar, decode_seq},
    fields::Fields,
};
use byteorder::{LittleEndian, ReadBytesExt};

use crate::{query::cigar::Cigar, query::cigar::Op, U32_SIZE};

#[cfg(not(feature = "python-ffi"))]
#[derive(Debug)]
/// Represents a GBAM record in which some fields may be omitted.
pub struct GbamRecord {
    /// Reference sequence ID
    pub refid: Option<i32>,
    /// 0-based leftmost coordinate
    pub pos: Option<i32>,
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
    pub cigar: Option<Cigar>,
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
    pub refid: Option<i32>,
    /// 0-based leftmost coordinate
    #[pyo3(get)]
    pub pos: Option<i32>,
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

/// cfg_attribute currently doesn't work with pyo3(get)
#[cfg(feature = "python-ffi")]
pub fn parse_cigar(bytes: &[u8]) -> String {
    decode_cigar(bytes)
}

/// This version is for Rust.
pub fn parse_cigar(bytes: &[u8]) -> Cigar {
    Cigar::new(
        bytes
            .chunks(U32_SIZE)
            .map(|mut slice| Op::new(slice.read_u32::<LittleEndian>().unwrap()))
            .collect(),
    )
}

// TODO :: ADD TEMPLATE LENGTHS TO GBAM RECORD
// TODO :: REMOVE CG TAG FROM ORIGINAL FILE
impl GbamRecord {
    pub(crate) fn parse_from_bytes(&mut self, field: &Fields, mut bytes: &[u8]) {
        match field {
            Fields::RefID => self.refid = Some(bytes.read_i32::<LittleEndian>().unwrap()),
            Fields::Pos => self.pos = Some(bytes.read_i32::<LittleEndian>().unwrap()),
            Fields::Mapq => self.mapq = Some(bytes[0].to_owned()),
            Fields::Bin => self.bin = Some(bytes.read_u16::<LittleEndian>().unwrap()),
            Fields::Flags => self.flag = Some(bytes.read_u16::<LittleEndian>().unwrap()),
            Fields::NextRefID => self.next_ref_id = Some(bytes.read_i32::<LittleEndian>().unwrap()),
            Fields::NextPos => self.next_pos = Some(bytes.read_i32::<LittleEndian>().unwrap()),
            Fields::TemplateLength => self.tlen = Some(bytes.read_i32::<LittleEndian>().unwrap()),
            Fields::ReadName => self.read_name = Some(bytes.to_vec()),
            Fields::RawCigar => self.cigar = Some(parse_cigar(bytes)),
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
