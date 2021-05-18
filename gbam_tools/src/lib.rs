//! This crate contains tools for interacting with GBAM file format.
#![deny(missing_docs)]

use self::Fields::*;

use std::mem;

use std::slice::Iter;

use serde::{Deserialize, Serialize};

/// BAM to GBAM converter
pub mod bam_to_gbam;
/// Meta information for GBAM file
pub mod meta;
/// GBAM reader
pub mod reader;
/// BAM record module
pub mod record;
/// Module responsible for tags parsing
mod tags;
/// GBAM writer
pub mod writer;

// use self::writer::Writer;
pub use self::reader::{ParsingTemplate, Reader};
use self::writer::Writer;
pub use crate::bam_to_gbam::bam_to_gbam;
pub use meta::Codecs;
pub use record::{decode_cigar, decode_seq, RawRecord};
use tags::get_tag;

const U64_SIZE: usize = mem::size_of::<u64>();
const U32_SIZE: usize = mem::size_of::<u32>();
const U16_SIZE: usize = mem::size_of::<u16>();
const U8_SIZE: usize = mem::size_of::<u8>();
const MEGA_BYTE_SIZE: usize = 1_048_576;

const SIZE_LIMIT: usize = 16 * MEGA_BYTE_SIZE;
static GBAM_MAGIC: &[u8] = b"geeBAM10";

const FIELDS_NUM: usize = 18;
/// Fields which contain data (not index fields).
const DATA_FIELDS_NUM: usize = 13;
/// Types of fields contained in BAM file.
#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[allow(missing_docs)]
pub enum Fields {
    /// Data fields
    #[allow(clippy::upper_case_acronyms)]
    RefID,
    Pos,
    Mapq,
    Bin,
    Flags,
    #[allow(clippy::upper_case_acronyms)]
    NextRefID,
    NextPos,
    TemplateLength,
    ReadName,
    RawCigar,
    RawSequence,
    RawQual,
    RawTags,
    /// Index fields
    LName,
    NCigar,
    SequenceLength,
    RawTagsLen, // Not in BAM spec, needed for index GBAM file
    RawSeqLen,  // Not in BAM spec, needed for index GBAM file
}

impl Fields {
    /// Returns iterator over enum fields
    pub fn iterator() -> Iter<'static, Fields> {
        static FIELDS: [Fields; FIELDS_NUM] = [
            // Data fields
            RefID,
            Pos,
            Mapq,
            Bin,
            Flags,
            NextRefID,
            NextPos,
            TemplateLength,
            ReadName,
            RawCigar,
            RawSequence,
            RawQual,
            RawTags,
            // Index fields
            LName,
            NCigar,
            SequenceLength,
            RawSeqLen,
            RawTagsLen,
        ];
        FIELDS.iter()
    }
}

/// Fields holding index are not data fields
pub(crate) fn is_data_field(field: &Fields) -> bool {
    matches!(
        field,
        RefID
            | Pos
            | Mapq
            | Bin
            | Flags
            | NextRefID
            | NextPos
            | TemplateLength
            | ReadName
            | RawCigar
            | RawSequence
            | RawQual
            | RawTags
    )
}

/// Variable sized fields have fixed size value NONE
/// This is **not** BAM sizes!
pub(crate) fn field_item_size(field: &Fields) -> Option<usize> {
    match field {
        RefID => Some(U32_SIZE),
        Pos => Some(U32_SIZE),
        LName => Some(U32_SIZE),
        Mapq => Some(U8_SIZE),
        Bin => Some(U16_SIZE),
        NCigar => Some(U32_SIZE),
        Flags => Some(U16_SIZE),
        SequenceLength => Some(U32_SIZE),
        NextRefID => Some(U32_SIZE),
        NextPos => Some(U32_SIZE),
        TemplateLength => Some(U32_SIZE),
        RawTagsLen => Some(U32_SIZE),
        RawSeqLen => Some(U32_SIZE),
        ReadName => None,
        RawCigar => None,
        RawSequence => None,
        RawQual => None,
        RawTags => None,
    }
}

impl std::fmt::Display for Fields {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        std::fmt::Debug::fmt(self, f)
    }
}

/// Type of Field.
#[derive(Debug)]
#[allow(missing_docs)]
pub(crate) enum FieldType {
    VariableSized,
    FixedSized,
}

// May be useful:
pub(crate) fn field_type(field: &Fields) -> FieldType {
    match field {
        // Fixed size fields
        Fields::RefID
        | Fields::Pos
        | Fields::LName
        | Fields::Mapq
        | Fields::Bin
        | Fields::NCigar
        | Fields::Flags
        | Fields::SequenceLength
        | Fields::NextRefID
        | Fields::NextPos
        | Fields::TemplateLength
        | Fields::RawSeqLen
        | Fields::RawTagsLen => FieldType::FixedSized,
        // Variable size fields
        Fields::ReadName
        | Fields::RawCigar
        | Fields::RawSequence
        | Fields::RawQual
        | Fields::RawTags => FieldType::VariableSized,
    }
}

/// Returns enum name of index field for particular variable sized field
fn var_size_field_to_index(field: &Fields) -> Fields {
    match field {
        Fields::ReadName => Fields::LName,
        Fields::RawQual => Fields::SequenceLength,
        Fields::RawSequence => Fields::RawSeqLen,
        Fields::RawTags => Fields::RawTagsLen,
        Fields::RawCigar => Fields::NCigar,
        _ => panic!("Unreachable"),
    }
}

use pyo3::prelude::*;
use pyo3::wrap_pyfunction;

/// Workaround, since it seems wrap_pyfunction cant access another module namespace.
#[pyfunction]
pub fn bam_to_gbam_python(in_path: String, out_path: String, codec_str: String) {
    let codec = match &codec_str[..] {
        "gzip" => Codecs::Gzip,
        "lz4" => Codecs::Lz4,
        _ => panic!("Codec <{}> is not supported.", codec_str),
    };
    bam_to_gbam(in_path, out_path, codec);
}

#[pymodule]
fn gbam_tools(_: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(bam_to_gbam_python, m)?)
        .unwrap();
    m.add_class::<reader::Reader>()?;
    m.add_class::<reader::ParsingTemplate>()?;
    m.add_class::<reader::GbamRecord>()?;
    Ok(())
}
