// https://github.com/pezmaster31/bamtools/blob/2391b1a1275816ad89c624586fa02b1a621924f5/src/api/internal/bam/BamReader_p.cpp

use crate::{U16_SIZE, U32_SIZE, U8_SIZE};
use byteorder::{LittleEndian, ReadBytesExt};

#[derive(Debug)]
pub(crate) enum TagType {
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
        _ => panic!("There is no tag type <{}>!\n", *c as char),
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
fn get_tag_data(data: &[u8]) -> (&[u8], usize, TagType) {
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
                // tag type + item type + byte count + len_in_bytes
                U8_SIZE + U8_SIZE + U32_SIZE + len_in_bytes,
                tag_type,
            )
        }
        TagType::Z | TagType::H => {
            while data[idx] != 0 {
                idx += 1;
            }
            (&data[U8_SIZE..idx], idx + 1, tag_type)
        }
        _ => {
            let item_size = tag_size(&tag_type).unwrap();
            (&data[idx..idx + item_size], U8_SIZE + item_size, tag_type)
        }
    }
}

pub(crate) fn get_tag<'a>(data: &'a [u8], tag: &[u8; 2]) -> Option<(&'a [u8], TagType)> {
    let mut idx: usize = 0;
    while idx < data.len() {
        let tag_value = get_tag_data(&data[idx + U16_SIZE..]);
        if &data[idx..idx + U16_SIZE] == tag {
            return Some((tag_value.0, tag_value.2));
        }
        // Skip tag
        idx += U16_SIZE;
        let tag_data_len = tag_value.1;
        // Skip tag value
        idx += tag_data_len;
    }
    None
}

// Returns value of HI tag.
// The field type is i so it's assumed it will fit in i32.
pub fn get_hit_count(data: &[u8]) -> Option<i32> {
    if let Some((mut tag, tag_type)) = get_tag(data, &[b'H', b'I']) {
        let val = match tag_type {
            TagType::A => tag.read_u8().unwrap() as i32,
            TagType::c => tag.read_i8().unwrap() as i32,
            TagType::C => tag.read_u8().unwrap() as i32,
            TagType::s => tag.read_i16::<LittleEndian>().unwrap() as i32,
            TagType::S => tag.read_u16::<LittleEndian>().unwrap() as i32,
            TagType::i => tag.read_i32::<LittleEndian>().unwrap(),
            TagType::I => tag.read_u32::<LittleEndian>().unwrap() as i32,
            _ => panic!("The tag type {:?} can't contain hit count value.", tag_type),
        };
        return Some(val);
    }
    None
}
