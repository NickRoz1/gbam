mod block;
pub mod gz;
mod util;
mod virtual_position;

// Module responsible for sorting.
pub mod sorting {
    mod comparators;
    // Contains flag parsing mechanism from noodles crate to extract flag data
    // from bamrawrecord.
    mod flags;
    pub mod sort;
}
mod reader;

pub mod record {
    /// BAM (raw) record module
    pub mod bamrawrecord;
    /// This module contains definition of Fields enum which is used to query
    /// BAM Raw Record fields.
    pub mod fields;
    /// Module responsible for tags parsing
    mod tags;
}

use block::Block;
pub use reader::parse_reference_sequences;
pub use reader::Reader;
use std::mem;
use virtual_position::VirtualPosition;

pub const MEGA_BYTE_SIZE: usize = 1024 * 1024;
#[allow(dead_code)]
pub(crate) const GIGA_BYTE_SIZE: usize = 1024 * MEGA_BYTE_SIZE;
#[allow(dead_code)]
const U64_SIZE: usize = mem::size_of::<u64>();
const U32_SIZE: usize = mem::size_of::<u32>();
const U16_SIZE: usize = mem::size_of::<u16>();
const U8_SIZE: usize = mem::size_of::<u8>();
const MAGIC_NUMBER: &[u8] = b"BAM\x01";
