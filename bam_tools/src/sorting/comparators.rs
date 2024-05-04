use super::flags;
use super::sort::SortBy;
use crate::record::bamrawrecord::BAMRawRecord;
use crate::record::fields::Fields;
use byteorder::{LittleEndian, ReadBytesExt};
use std::borrow::Cow;
use std::cmp::Ordering;
use std::str::from_utf8_unchecked;

/// This enum hold variants which hold keys for each of the available sort
/// variations. It is done to avoid repeated calling of key extracting functions
/// (like get hit count tag), sort of a preliminary caching (proper caching is
/// not trivial to implement in this case).
#[derive(Debug)]
pub(crate) enum KeyTuple<'a> {
    // Holds name and original record (to ensure that name is valid throughout time)
    Name(Cow<'a, str>),
    // Holds name, hit count, flag and original record (to ensure that name is valid throughout time)
    NameAndMatchMates(Cow<'a, str>, Option<i32>, u16),
    // Holds RefID, POS and is_reversed_strand
    CoordinatesAndStrand(i32, i32, bool),
}

// pub(crate) fn extract_key_alloc_on_arena<'a>(
//     rec: &BAMRawRecord,
//     sort_by: &SortBy,
//     read_names_arena: &'a Arena<String>,
// ) -> KeyTuple<'a> {
//     let slice = rec.get_bytes(&Fields::ReadName);
//     let str = unsafe { from_utf8_unchecked(slice) };
//     let name = read_names_arena.alloc(str.to_owned());
//     create_key_tuple(Cow::Borrowed(name), rec, sort_by)
// }

pub(crate) fn extract_key<'a>(rec: &BAMRawRecord, buf: &'a [u8], sort_by: &SortBy) -> KeyTuple<'a> {
    let slice = rec.get_range(&Fields::ReadName);
    let name = unsafe { from_utf8_unchecked(&buf[slice]) };
    create_key_tuple(name, rec, sort_by)
}

pub(crate) fn create_key_tuple<'a>(
    name: &'a str,
    rec: &BAMRawRecord,
    sort_by: &SortBy,
) -> KeyTuple<'a> {
    match *sort_by {
        SortBy::Name => KeyTuple::Name(Cow::Borrowed(name)),
        SortBy::NameAndMatchMates => {
            KeyTuple::NameAndMatchMates(Cow::Borrowed(name), rec.get_hit_count(), get_flag_val(rec))
        }
        SortBy::CoordinatesAndStrand => {
            let mut id = get_ref_id(rec);
            if id == -1 {
                id = i32::MAX;
            }
            KeyTuple::CoordinatesAndStrand(id, get_pos(rec), is_reverse_strand(rec))
        }
    }
}

pub(crate) fn compare_read_names(lhs: &KeyTuple, rhs: &KeyTuple) -> Ordering {
    let left = get_name(lhs);
    let right = get_name(rhs);
    left.cmp(right)
}

fn get_name<'a>(key_tuple: &'a KeyTuple) -> &'a str {
    match &key_tuple {
        KeyTuple::Name(str) => str,
        KeyTuple::NameAndMatchMates(str, _, _) => str,
        _ => panic!("Variant {:?} does not have Name field.", key_tuple),
    }
}

/// Comparison function for 'queryname' sorting order setting mates
/// of the same alignments adjacent with the first mate coming before
/// the second mate
pub(crate) fn compare_read_names_and_mates(lhs: &KeyTuple, rhs: &KeyTuple) -> Ordering {
    let ordering = compare_read_names(lhs, rhs);
    if ordering == Ordering::Equal {
        let tag_lhs = extract_hit_count(lhs);
        let tag_rhs = extract_hit_count(rhs);

        if tag_lhs == tag_rhs {
            // Source https://github.com/biod/sambamba/blob/c795656721b3608ffe7765b6ab98502426d14131/BioD/bio/std/hts/bam/read.d#L1603
            let flag_lhs = extract_flag(lhs);
            let flag_rhs = extract_flag(rhs);
            return flag_lhs.cmp(flag_rhs);
        }
        return tag_lhs.unwrap().cmp(&tag_rhs.unwrap());
    }
    ordering
}

fn extract_hit_count<'a>(key_tuple: &'a KeyTuple) -> &'a Option<i32> {
    if let KeyTuple::NameAndMatchMates(_, hit_count, _) = key_tuple {
        hit_count
    } else {
        panic!("Variant {:?} does not have Hit Count field.", key_tuple);
    }
}

fn extract_flag<'a>(key_tuple: &'a KeyTuple) -> &'a u16 {
    if let KeyTuple::NameAndMatchMates(_, _, flag) = key_tuple {
        flag
    } else {
        panic!("Variant {:?} does not have Flag field.", key_tuple);
    }
}

fn get_flag_val(rec: &BAMRawRecord) -> u16 {
    rec.get_bytes(&Fields::Flags)
        .read_u16::<LittleEndian>()
        .unwrap()
}

pub(crate) fn compare_coordinates_and_strand(left: &KeyTuple, right: &KeyTuple) -> Ordering {
    if let (KeyTuple::CoordinatesAndStrand(_, _, _), KeyTuple::CoordinatesAndStrand(_, _, _)) =
        (left, right)
    {
    } else {
        panic!(
            "CoordinatesAndStrand variant is expected. Got {:?} and {:?}",
            left, right
        );
    }
    let refid_left = extract_ref_id(left);
    let refid_right = extract_ref_id(right);
    if refid_left != refid_right {
        return refid_left.cmp(refid_right);
    }
    let pos_left = extract_pos(left);
    let pos_right = extract_pos(right);
    if pos_left != pos_right {
        return pos_left.cmp(pos_right);
    }
    let is_reverse_strand_left = extract_is_reversed(left);
    let is_reverse_strand_right = extract_is_reversed(right);

    if !is_reverse_strand_left && *is_reverse_strand_right {
        Ordering::Less
    } else if *is_reverse_strand_left && !is_reverse_strand_right {
        Ordering::Greater
    } else {
        Ordering::Equal
    }
}

fn extract_ref_id<'a>(key_tuple: &'a KeyTuple) -> &'a i32 {
    if let KeyTuple::CoordinatesAndStrand(ref_id, _, _) = key_tuple {
        ref_id
    } else {
        panic!("Variant {:?} does not have RefId field.", key_tuple);
    }
}

fn extract_pos<'a>(key_tuple: &'a KeyTuple) -> &'a i32 {
    if let KeyTuple::CoordinatesAndStrand(_, pos, _) = key_tuple {
        pos
    } else {
        panic!("Variant {:?} does not have Pos field.", key_tuple);
    }
}

fn extract_is_reversed<'a>(key_tuple: &'a KeyTuple) -> &'a bool {
    if let KeyTuple::CoordinatesAndStrand(_, _, is_reversed) = key_tuple {
        is_reversed
    } else {
        panic!("Variant {:?} does not have Is Reversed field.", key_tuple);
    }
}

fn is_reverse_strand(rec: &BAMRawRecord) -> bool {
    let flags = rec
        .get_bytes(&Fields::Flags)
        .read_u16::<LittleEndian>()
        .unwrap();
    let bit_field = flags::Flags::from_bits(flags).unwrap();
    bit_field.is_reverse_complemented()
}

fn get_ref_id(rec: &BAMRawRecord) -> i32 {
    rec.get_bytes(&Fields::RefID)
        .read_i32::<LittleEndian>()
        .unwrap()
}

fn get_pos(rec: &BAMRawRecord) -> i32 {
    rec.get_bytes(&Fields::Pos)
        .read_i32::<LittleEndian>()
        .unwrap()
}
