use std::convert::TryInto;
use std::ops::{Range, RangeInclusive};
use std::{cmp::Ordering, collections::HashSet, convert::TryFrom, time::Instant};

use bam_tools::record::fields::Fields;

/// This module provides function for fast querying of read depth.
use crate::meta::BlockMeta;
use crate::reader::{reader::generate_block_treemap, reader::Reader, record::GbamRecord};

// use crate::reader::record::copying;

use byteorder::{LittleEndian, ReadBytesExt};

// Approach as in https://github.com/brentp/mosdepth
pub fn get_region_depth(region: (&str, i32, i32)) {
    if !reader
        .parsing_template
        .check_if_active(&[Fields::RefID, Fields::Pos, Fields::RawCigar])
    {
        panic!("The reader should have parsing template which includes REFID, POS and CIGAR.");
    }

    let meta_blocks = reader.file_meta.view_blocks(&Fields::RefID);

    let mut left_bound = 0;
    let mut right_bound = reader.amount;

    // If stats were collected for this file it's possible to narrow binary search.
    if meta_blocks.first().unwrap().max_value.is_some() {
        // println!("yeah that is stat");
        // First block which contains this refid.
        let first_block = find_leftmost_block(chr_num, &meta_blocks);

        // There is no such refid inside.
        if first_block == meta_blocks.len() {
            return None;
        }

        left_bound = meta_blocks[..first_block]
            .iter()
            .fold(0, |x, meta| x + meta.numitems) as usize;
        let block_size = meta_blocks[first_block].numitems;
        right_bound = left_bound + block_size as usize;
    }

    // Find number of record when records with this ref_id begin.
    // We don't need to fetch anything except refid to find first record with requested refid, so disable other fields fetching.
    reader.fetch_only(&[Fields::RefID]);
    // println!("Look - the template: {:?}", reader.parsing_template);
    // println!("Left: {}, Right: {}", left_bound, right_bound);
    let res = find_refid(reader, chr_num, left_bound..right_bound);
    reader.restore_template();

    if let Some(region_start) = res {
        Some(calc_depth_at_pos(reader, region_start, pos, chr_num))
    } else {
        None
    }
}

/// Get depth at position. The file should be sorted!
pub fn get_depth(reader: &mut Reader, chr: String, pos: i32) -> Option<usize> {
    if !reader
        .parsing_template
        .check_if_active(&[Fields::RefID, Fields::Pos, Fields::RawCigar])
    {
        panic!("The reader should have parsing template which includes REFID, POS and CIGAR.");
    }

    // https://github.com/biod/sambamba/blob/3eff9a2d8bb3097b92c72752be3c6b42dd1c59b7/BioD/bio/std/hts/bam/read.d#L902
    let chr_num = match reader.file_meta.get_ref_id(chr) {
        Some(num) => num,
        None => return None,
    };

    // // Used to quickly determine what block record belongs to.
    // let blocks = generate_block_treemap(&reader.file_meta, &Fields::RawCigar);

    let meta_blocks = reader.file_meta.view_blocks(&Fields::RefID);

    let mut left_bound = 0;
    let mut right_bound = reader.amount;

    // If stats were collected for this file it's possible to narrow binary search.
    if meta_blocks.first().unwrap().max_value.is_some() {
        // println!("yeah that is stat");
        // First block which contains this refid.
        let first_block = find_leftmost_block(chr_num, &meta_blocks);

        // There is no such refid inside.
        if first_block == meta_blocks.len() {
            return None;
        }

        left_bound = meta_blocks[..first_block]
            .iter()
            .fold(0, |x, meta| x + meta.numitems) as usize;
        let block_size = meta_blocks[first_block].numitems;
        right_bound = left_bound + block_size as usize;
    }

    // Find number of record when records with this ref_id begin.
    // We don't need to fetch anything except refid to find first record with requested refid, so disable other fields fetching.
    reader.fetch_only(&[Fields::RefID]);
    // println!("Look - the template: {:?}", reader.parsing_template);
    // println!("Left: {}, Right: {}", left_bound, right_bound);
    let res = find_refid(reader, chr_num, left_bound..right_bound);
    reader.restore_template();

    if let Some(region_start) = res {
        Some(calc_depth_at_pos(reader, region_start, pos, chr_num))
    } else {
        None
    }
}

fn find_leftmost_block(id: i32, block_metas: &Vec<BlockMeta>) -> usize {
    let mut left = 0;
    let mut right = block_metas.len();
    while left < right {
        let mid = (left + right) / 2;
        let max_val = (&block_metas[mid].max_value.as_ref().unwrap()[..])
            .read_i32::<LittleEndian>()
            .unwrap();
        match max_val.cmp(&id) {
            Ordering::Equal | Ordering::Greater => right = mid,
            Ordering::Less => left = mid + 1,
        }
    }
    left
}

// For each guess there may be I/O operation with decompression, so this method is not fast.
fn find_refid(reader: &mut Reader, chr_num: i32, range: Range<usize>) -> Option<usize> {
    let pred = |num: usize, buf: &mut GbamRecord| {
        reader.fill_record(num, buf);
        buf.refid.unwrap().cmp(&chr_num)
    };

    binary_search(range.start, range.end, pred)
}

/// Searches for the first record which satisfies predicate.
fn binary_search<F: FnMut(usize, &mut GbamRecord) -> Ordering>(
    mut left: usize,
    mut right: usize,
    mut cmp: F,
) -> Option<usize> {
    let mut buf = GbamRecord::default();
    let end = right;

    while left < right {
        let mid = (left + right) / 2;
        // println!("Mid - {:?}", mid);
        match cmp(mid, &mut buf) {
            Ordering::Less => left = mid + 1,
            Ordering::Equal | Ordering::Greater => right = mid,
        }
    }

    if left == end || cmp(left, &mut buf) != Ordering::Equal {
        return None;
    }
    Some(left)
}

pub fn calc_depth(
    reader: &mut Reader,
    mut record_num: usize,
    region: RangeInclusive<i32>,
    refid: i32,
) -> Vec<i32> {
    let mut buf = GbamRecord::default();
    let mut sweeping_line = vec![0; (region.end() - region.start() + 1) as usize];

    loop {
        if record_num >= reader.amount {
            break;
        }
        reader.fill_record(record_num, &mut buf);
        record_num += 1;

        let read_start = buf.pos.unwrap();
        let base_cov = buf.cigar.as_ref().unwrap().base_coverage();

        if read_start >= *region.end() || buf.refid.unwrap() > refid {
            break;
        }

        let read_end = read_start + i32::try_from(base_cov).unwrap();

        if read_start > *region.start() {
            sweeping_line[(read_start - region.start()) as usize] += 1;
        }
        if read_end > *region.start() && read_end < *region.end() {
            sweeping_line[(read_start - region.start()) as usize] -= 1;
        }
    }

    let mut depth_of_region = sweeping_line;
    let mut cur_depth = 0;
    let mut delta = 0;

    for ref mut sweep in depth_of_region {
        debug_assert!(cur_depth >= 0);
        delta = *sweep;
        *sweep = cur_depth;
        cur_depth += delta;
    }

    return depth_of_region;
}

// println!(
//     "REF_ID: {}\nPOS: {}\nbase_cov: {}\ncigar_len: {}\nsearch pos: {}\nEND: {}\nDEPTH: {}\n",
//     buf.refid.unwrap(),
//     read_pos,
//     base_cov,
//     buf.cigar.as_ref().unwrap().0.len(),
//     pos,
//     read_pos + i32::try_from(base_cov).unwrap(),
//     read_pos > pos
// );
