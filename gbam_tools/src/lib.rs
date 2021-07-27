//! This crate contains tools for interacting with GBAM file format.
#![allow(missing_docs)]

use std::mem;

pub mod bam {
    /// BAM to GBAM converter
    pub mod bam_to_gbam;
}

pub mod reader {
    mod column;
    /// GBAM reader
    pub mod reader;
}
/// Meta information for GBAM file
pub mod meta;
/// GBAM writer
pub mod writer;

mod compressor;

// use self::writer::Writer;
// pub use {ParsingTemplate, Reader};
use self::writer::Writer;
pub use bam::bam_to_gbam::{bam_sort_to_gbam, bam_to_gbam};
use compressor::{CompressTask, Compressor};
pub use meta::Codecs;

const U64_SIZE: usize = mem::size_of::<u64>();
const U32_SIZE: usize = mem::size_of::<u32>();
const U16_SIZE: usize = mem::size_of::<u16>();
const U8_SIZE: usize = mem::size_of::<u8>();
const MEGA_BYTE_SIZE: usize = 1_048_576;

/// 16777216 bytes
const SIZE_LIMIT: usize = 16 * MEGA_BYTE_SIZE;
static GBAM_MAGIC: &[u8] = b"geeBAM10";

#[cfg(feature = "python-ffi")]
mod ffi {
    use crate::reader;
    use crate::Codecs;
    use crate::{bam_sort_to_gbam, bam_to_gbam};

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
        m.add_class::<reader::Reader>()?;
        m.add_class::<reader::ParsingTemplate>()?;
        m.add_class::<reader::GbamRecord>()?;
        Ok(())
    }
}
