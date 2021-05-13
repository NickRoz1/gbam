// https://github.com/pezmaster31/bamtools/blob/2391b1a1275816ad89c624586fa02b1a621924f5/src/api/internal/bam/BamReader_p.cpp

use super::{U16_SIZE, U32_SIZE, U8_SIZE};
use byteorder::{LittleEndian, ReadBytesExt};

enum TagType {
    /// Char
    A,
    /// Byte array
    B,
    /// u8
    C,
    /// i8
    #[allow(non_camel_case_types)]
    c,
    /// float
    #[allow(non_camel_case_types)]
    f,
    /// Null-terminated HEX string
    H,
    /// u32
    I,
    /// i32
    #[allow(non_camel_case_types)]
    i,
    /// u16
    S,
    /// i16
    #[allow(non_camel_case_types)]
    s,
    /// Null-terminated char string
    Z,
}

fn get_tag_type(c: &u8) -> TagType {
    match *c as char {
        'A' => TagType::A,
        'B' => TagType::B,
        'C' => TagType::C,
        'c' => TagType::c,
        'f' => TagType::f,
        'H' => TagType::H,
        'i' => TagType::i,
        'I' => TagType::I,
        'S' => TagType::S,
        's' => TagType::s,
        'Z' => TagType::Z,
        _ => panic!("There is no tag type <{}>!\n", c),
    }
}

fn tag_size(tag: &TagType) -> Option<usize> {
    match *tag {
        TagType::C | TagType::c | TagType::A => Some(U8_SIZE),
        TagType::S | TagType::s => Some(U16_SIZE),
        TagType::I | TagType::i | TagType::f => Some(U32_SIZE),
        _ => None,
    }
}

/// Returns slice with tag data and amount it and related info occupies in
/// original buffer
fn get_tag_data(data: &[u8]) -> (&[u8], usize) {
    let mut idx = 0;
    let tag_type = get_tag_type(&data[idx]);
    idx += U8_SIZE;
    match tag_type {
        TagType::B => {
            let item_size = tag_size(&get_tag_type(&data[idx])).unwrap();
            idx += U8_SIZE;
            let mut len_bytes = &data[idx..idx + U32_SIZE];
            let len = len_bytes.read_u32::<LittleEndian>().unwrap() as usize;
            let len_in_bytes = len * item_size;
            idx += U32_SIZE;
            (
                &data[idx..idx + len_in_bytes],
                U8_SIZE + U8_SIZE + U32_SIZE + len_in_bytes,
            )
        }
        TagType::Z | TagType::H => {
            while data[idx] != 0 {
                idx += 1;
            }
            (&data[U8_SIZE..idx], idx)
        }
        _ => {
            let item_size = tag_size(&tag_type).unwrap();
            (&data[idx..idx + item_size], U8_SIZE + item_size)
        }
    }
}

pub(crate) fn get_tag<'a>(data: &'a [u8], tag: &[u8]) -> Option<&'a [u8]> {
    let mut idx: usize = 0;
    while idx < data.len() {
        if &data[idx..idx + U16_SIZE] == tag {
            return Some(get_tag_data(&data[idx + U16_SIZE..]).0);
        }
        // Skip tag
        idx += U16_SIZE;
        let tag_data_len = get_tag_data(&data[idx..]).1;
        // Skip tag value
        idx += tag_data_len;
    }
    None
}
