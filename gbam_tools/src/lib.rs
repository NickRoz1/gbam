//! This crate contains tools for interacting with GBAM file format.
#![deny(missing_docs)]

use self::Fields::*;
use byteorder::{ByteOrder, LittleEndian};
use std::mem;
use std::ops::{Deref, DerefMut};
use std::slice::Iter;

mod compression;
mod meta;
// mod reader;
/// BAM to GBAM converter
pub mod bam_to_gbam;
mod rowgroup;
/// Single threaded reader
mod single_thread;
/// Writer module
mod writer;

// use self::writer::Writer;
use self::single_thread::writer::Writer;
pub use crate::bam_to_gbam::bam_to_gbam;
use crate::compression::{Compression, COMPRESSION_ENUM_SIZE};
use crate::meta::{ColChunkMeta, RowGroupMeta};

static GBAM_MAGIC: &[u8] = b"GBAM-0.1.0\x01";

const u64_size: usize = mem::size_of::<u64>();
const u32_size: usize = mem::size_of::<u32>();
const u16_size: usize = mem::size_of::<u16>();
const u8_size: usize = mem::size_of::<u8>();
const mega_byte_size: usize = 1_048_576;

const FIELDS_NUM: usize = 17;
/// Types of fields contained in BAM file.
#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
#[allow(missing_docs)]
pub enum Fields {
    RefID,
    Pos,
    LName,
    Mapq,
    Bin,
    NCigar,
    Flags,
    SequenceLength,
    NextRefID,
    NextPos,
    TemplateLength,
    ReadName,
    RawCigar,
    RawSequence,
    RawQual,
    RawTags,
    RawTagsLen, // Not in BAM spec, needed for index GBAM file
    RawSeqLen,  // Not in BAM spec, needed for index GBAM file
}

impl Fields {
    /// Returns iterator over enum fields
    pub fn iterator() -> Iter<'static, Fields> {
        static FIELDS: [Fields; FIELDS_NUM] = [
            RefID,
            Pos,
            LName,
            Mapq,
            Bin,
            NCigar,
            Flags,
            SequenceLength,
            NextRefID,
            NextPos,
            TemplateLength,
            ReadName,
            RawCigar,
            RawSequence,
            RawQual,
            RawTags,
            RawTagsLen,
        ];
        FIELDS.iter()
    }
}

// match field {
//     // Fixed size fields
//     Fields::RefID
//     | Fields::Pos
//     | Fields::LName
//     | Fields::Mapq
//     | Fields::Bin
//     | Fields::NCigar
//     | Fields::Flags
//     | Fields::SequenceLength
//     | Fields::NextRefID
//     | Fields::NextPos
//     | Fields::TemplateLength
//     | Fields::RawTagsLen => {}
//     // Variable size fields
//     Fields::ReadName
//     | Fields::RawCigar
//     | Fields::RawSequence
//     | Fields::RawQual
//     | Fields::RawTags => {
//     }
// }

/// Provides convenient access to record bytes
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct RawRecord(Vec<u8>);

impl RawRecord {
    /// Change size of record
    pub fn resize(&mut self, new_len: usize) {
        self.0.resize(new_len, Default::default());
    }

    fn get_slice(&self, offset: usize, len: usize) -> &[u8] {
        &self.0[offset..offset + len]
    }

    fn l_read_name(&self) -> u8 {
        self.get_bytes(&Fields::LName)[0]
    }

    fn n_cigar_op(&self) -> u16 {
        LittleEndian::read_u16(self.get_bytes(&Fields::NCigar))
    }

    fn l_seq(&self) -> u32 {
        LittleEndian::read_u32(self.get_bytes(&Fields::SequenceLength))
    }

    /// Values of fields containg length of other fields
    pub fn get_len_val(&self, field: &Fields) -> usize {
        match field {
            LName => self.l_read_name() as usize,
            SequenceLength => self.l_seq() as usize,
            NCigar => self.n_cigar_op() as usize,
            RawTagsLen => self.0.len() - self.get_offset(&Fields::RawTags),
            _ => panic!("This field is not supported: {} \n", *field as usize),
        }
    }

    /// Calculates actual size of variable length field in bytes.
    pub fn get_var_field_len(&self, field: &Fields) -> usize {
        match field {
            ReadName => self.l_read_name() as usize,
            RawCigar => u32_size * self.n_cigar_op() as usize,
            RawSequence => ((self.l_seq() + 1) / 2) as usize,
            RawQual => self.l_seq() as usize,
            RawTags => self.0.len() - self.get_offset(&Fields::RawTags),
            _ => panic!("This field is not supported: {} \n", *field as usize),
        }
    }

    fn get_offset(&self, field: &Fields) -> usize {
        match field {
            ReadName => 32,
            RawCigar => {
                self.get_offset(&Fields::ReadName) + self.get_var_field_len(&Fields::ReadName)
            }
            RawSequence => {
                self.get_offset(&Fields::RawCigar) + self.get_var_field_len(&Fields::RawCigar)
            }
            RawQual => {
                self.get_offset(&Fields::RawSequence) + self.get_var_field_len(&Fields::RawSequence)
            }
            RawTags => self.get_offset(&Fields::RawQual) + self.get_var_field_len(&Fields::RawQual),
            _ => panic!("This field is not supported: {} \n", *field as usize),
        }
    }

    fn get_bytes(&self, field: &Fields) -> &[u8] {
        let get_cigar_offset = || -> usize { (32 + self.l_read_name()) as usize };
        let get_seq_offset =
            || -> usize { get_cigar_offset() + u32_size * self.n_cigar_op() as usize };
        let get_qual_offset = || -> usize { get_seq_offset() + ((self.l_seq() + 1) / 2) as usize };
        let get_tags_offset = || -> usize { get_qual_offset() + self.l_seq() as usize };
        match field {
            RefID => self.get_slice(0, u32_size),
            Pos => self.get_slice(4, u32_size),
            LName => self.get_slice(8, u8_size),
            Mapq => self.get_slice(9, u8_size),
            Bin => self.get_slice(10, u16_size),
            NCigar => self.get_slice(12, u16_size),
            Flags => self.get_slice(14, u16_size),
            SequenceLength => self.get_slice(16, u32_size),
            NextRefID => self.get_slice(20, u32_size),
            NextPos => self.get_slice(24, u32_size),
            TemplateLength => self.get_slice(28, u32_size),
            ReadName => self.get_slice(32, self.get_var_field_len(field)),
            RawCigar => self.get_slice(get_cigar_offset(), self.get_var_field_len(field)),
            RawSequence => self.get_slice(get_seq_offset(), self.get_var_field_len(field)),
            RawQual => self.get_slice(get_qual_offset(), self.l_seq() as usize),
            RawTags => self.get_slice(get_tags_offset(), self.0.len() - get_tags_offset()),
            _ => panic!("This field is not supported: {} \n", *field as usize),
        }
    }
}

impl From<Vec<u8>> for RawRecord {
    fn from(bytes: Vec<u8>) -> Self {
        Self(bytes)
    }
}

// Source: https://github.com/zaeleus/noodles/blob/316ec6f42960e4540bb2acc45b5653fb00b9970c/noodles-bam/src/record.rs#L324
impl Default for RawRecord {
    fn default() -> Self {
        Self::from(vec![
            0xff, 0xff, 0xff, 0xff, // ref_id = -1
            0xff, 0xff, 0xff, 0xff, // pos = -1
            0x02, // l_read_name = 2
            0xff, // mapq = 255
            0x48, 0x12, // bin = 4680
            0x00, 0x00, // n_cigar_op = 0
            0x04, 0x00, // flag = 4
            0x00, 0x00, 0x00, 0x00, // l_seq = 0
            0xff, 0xff, 0xff, 0xff, // next_ref_id = -1
            0xff, 0xff, 0xff, 0xff, // next_pos = -1
            0x00, 0x00, 0x00, 0x00, // tlen = 0
            0x2a, 0x00, // read_name = "*\x00"
        ])
    }
}

impl Deref for RawRecord {
    type Target = [u8];

    fn deref(&self) -> &[u8] {
        &self.0
    }
}

impl DerefMut for RawRecord {
    fn deref_mut(&mut self) -> &mut [u8] {
        &mut self.0
    }
}

use pyo3::prelude::*;
use pyo3::{wrap_pyfunction, wrap_pymodule};

/// Workaround, since it seems wrap_pyfunction cant access another module namespace.
#[pyfunction]
pub fn bam_to_gbam_python(in_path: String, out_path: String) {
    bam_to_gbam(in_path, out_path);
}

#[pymodule]
fn gbam_tools(_: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(bam_to_gbam_python, m)?)
        .unwrap();
    m.add_class::<single_thread::reader::Reader>()?;
    m.add_class::<single_thread::reader::ParsingTemplate>()?;
    m.add_class::<single_thread::reader::GbamRecord>()?;
    Ok(())
}
