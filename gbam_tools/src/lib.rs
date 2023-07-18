//! This crate contains tools for interacting with GBAM file format.
#![allow(missing_docs)]

use std::mem;

pub mod bam {
    /// BAM to GBAM converter
    pub mod bam_to_gbam;
    /// GBAM to BAM converter
    pub mod gbam_to_bam;
}

pub mod utils {
    /// BED reader
    pub mod bed;
}

pub mod reader {
    mod column;
    pub mod parse_tmplt;
    /// GBAM reader
    #[allow(clippy::module_inception)]
    pub mod reader;
    pub mod record;
    pub mod records;
}

#[cfg(not(feature = "python-ffi"))]
pub mod query {
    pub mod cigar;
    pub mod depth;
    pub mod flagstat;
    pub mod int2str;
}

/// Manages parallel compression
mod compressor;
/// Meta information for GBAM file
pub mod meta;
/// Manages stats collection
mod stats;
/// GBAM writer
pub mod writer;

// use self::writer::Writer;
// pub use {ParsingTemplate, Reader};
use self::writer::Writer;
pub use bam::bam_to_gbam::{bam_sort_to_gbam, bam_to_gbam};
pub use meta::Codecs;
pub use bam_tools::record::fields::Fields;


const U32_SIZE: usize = mem::size_of::<u32>();
const MEGA_BYTE_SIZE: usize = 1_048_576;

/// 16777216 bytes
const SIZE_LIMIT: usize = 8 * MEGA_BYTE_SIZE;
static GBAM_MAGIC: &[u8] = b"geeBAM10";

#[cfg(feature = "python-ffi")]
mod ffi {
    use crate::{
        Codecs, {bam_sort_to_gbam, bam_to_gbam},
    };

    use crate::reader::{parse_tmplt, record, records};

    use pyo3::prelude::*;
    use pyo3::wrap_pyfunction;

    /// Workaround, since it seems wrap_pyfunction cant access another module namespace.
    #[pyfunction]
    pub fn bam_to_gbam_python(in_path: String, out_path: String, codec_str: String, sort: bool) {
        let codec = match &codec_str[..] {
            "gzip" => Codecs::Gzip,
            "lz4" => Codecs::Lz4,
            _ => panic!("Codec <{}> is not supported.", codec_str),
        };
        if sort {
            bam_sort_to_gbam(&in_path, &out_path, codec);
        } else {
            bam_to_gbam(&in_path, &out_path, codec);
        }
    }

    #[pyfunction]
    pub fn test() {}

    #[pymodule]
    fn gbam_tools(_: Python, m: &PyModule) -> PyResult<()> {
        m.add_function(wrap_pyfunction!(test, m)?).unwrap();
        m.add_function(wrap_pyfunction!(bam_to_gbam_python, m)?)
            .unwrap();
        m.add_class::<records::PyRecords>()?;
        m.add_class::<parse_tmplt::ParsingTemplate>()?;
        m.add_class::<record::GbamRecord>()?;
        Ok(())
    }
}
