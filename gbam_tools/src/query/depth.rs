use std::convert::TryInto;
use std::ops::{Range, RangeInclusive};
use std::{cmp::Ordering, collections::HashSet, convert::TryFrom, time::Instant};

use bam_tools::record::fields::Fields;

/// This module provides function for fast querying of read depth.
use crate::meta::{BlockMeta, Limits};
use crate::reader::{reader::generate_block_treemap, reader::Reader, record::GbamRecord};

// use crate::reader::record::copying;

use byteorder::{LittleEndian, ReadBytesExt};

/// To avoid wasting buffer when get_region_depth returns None.
pub enum DepthStatus {
    Success(Vec<i32>),
    None(Vec<i32>),
}

// TODO: Merge bed regions with tolerance to form super regions. Then do
// sweepline for each of the super regions. Print the depth results for each of
// the regions by taking chunks from the super regions or doing in process while
// calculating to avoid double iterating the super region sweep line array.

fn parse_bytes(bytes: &Vec<u8>) -> i32 {
    (&bytes[..]).read_i32::<LittleEndian>().unwrap()
}

impl Limits<i32> for BlockMeta {
    fn get_max(&self) -> Option<i32> {
        self.max_value.and_then(|bytes| Some(parse_bytes(&bytes)))
    }
    fn get_min(&self) -> Option<i32> {
        self.min_value.and_then(|bytes| Some(parse_bytes(&bytes)))
    }
}

// Approach as in https://github.com/brentp/mosdepth
/// Get depth at position. The file should be sorted!
pub fn get_region_depths(
    reader: &mut Reader,
    regions: &Vec<&(String, i32, i32)>,
    buffer: Vec<i32>,
) -> DepthStatus {
    if !reader
        .parsing_template
        .check_if_active(&[Fields::RefID, Fields::Pos, Fields::RawCigar])
    {
        panic!("The reader should have parsing template which includes REFID, POS and CIGAR.");
    }

    // dbg!(&region.0);
    // https://github.com/biod/sambamba/blob/3eff9a2d8bb3097b92c72752be3c6b42dd1c59b7/BioD/bio/std/hts/bam/read.d#L902
    let chr_nums = regions
        .into_iter()
        .map(|region| reader.file_meta.get_ref_id(&region.0))
        .collect::<Vec<Option<i32>>>();

    let meta_blocks = reader.file_meta.view_blocks(&Fields::RefID);

    let mut left_bound = 0;
    let mut right_bound = reader.amount;

    let mut chrs_bounds: Vec<Option<(usize, usize)>>;
    chrs_bounds.resize(chr_nums.len(), None);

    // If stats were collected for this file it's possible to narrow binary search.
    if meta_blocks.first().unwrap().max_value.is_some() {
        let meta_iter = meta_blocks.iter();
        let mut cur_offset: usize = 0;
        for (i, cur_chr) in chr_nums.iter().flatten().enumerate() {
            // Iterate thorugh blocks, until one with max value bigger than cur_chr is found.
            while let Some(meta_block) = meta_iter.next() {
                if &meta_block.get_max().unwrap() >= cur_chr {
                    if &meta_block.get_min().unwrap() > cur_chr {
                        chrs_bounds[i] = None;
                    } else {
                        // This block may contain the REFID. Later on we will search this block to determine this ultimately.
                        chrs_bounds[i] =
                            Some((cur_offset, cur_offset + meta_block.numitems as usize));
                    }
                    break;
                }
                cur_offset += meta_block.numitems as usize;
            }
        }
    }

    for bound in chrs_bounds.iter().flatten() {}

    // Find number of record when records with this ref_id begin.
    // We don't need to fetch anything except refid to find first record with requested refid, so disable other fields fetching.
    reader.fetch_only(&[Fields::RefID]);
    let res = find_refid(reader, chr_num, left_bound..right_bound);
    reader.restore_template();

    if let Some(region_start) = res {
        calc_depth(reader, region_start, chr_num, region.1..=region.2, buffer)
    } else {
        DepthStatus::None(buffer)
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
    // dbg!("Looking for", chr_num);
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
        // dbg!(mid, &buf.refid);
        match cmp(mid, &mut buf) {
            Ordering::Less => left = mid + 1,
            Ordering::Equal | Ordering::Greater => right = mid,
        }
    }
    // dbg!("Finished search");
    if left == end || cmp(left, &mut buf) != Ordering::Equal {
        return None;
    }
    Some(left)
}

pub fn calc_depth(
    reader: &mut Reader,
    mut record_num: usize,
    refid: i32,
    region: RangeInclusive<i32>,
    buffer: Vec<i32>,
) -> DepthStatus {
    let mut buf = GbamRecord::default();

    let mut sweeping_line = buffer;
    sweeping_line.clear();
    sweeping_line.resize((region.end() - region.start() + 1) as usize, 0);
    let mut cur_depth = 0;

    loop {
        if record_num >= reader.amount {
            break;
        }
        reader.fill_record(record_num, &mut buf);
        record_num += 1;

        let read_start = buf.pos.unwrap();
        let base_cov = buf.cigar.as_ref().unwrap().base_coverage();

        let read_end = read_start + i32::try_from(base_cov).unwrap();

        // println!("That is test! {} end {}", read_start, read_end);

        if read_start >= *region.end() || buf.refid.unwrap() > refid {
            break;
        }

        if read_start >= *region.start() {
            sweeping_line[(read_start - region.start()) as usize] += 1;
        }
        if read_end > *region.start() && read_end <= *region.end() {
            if read_start < *region.start() {
                cur_depth += 1;
            }
            sweeping_line[(read_end - region.start()) as usize] -= 1;
        }
    }

    let mut depth_of_region = sweeping_line;
    let mut delta;

    for sweep in depth_of_region.iter_mut() {
        debug_assert!(cur_depth >= 0);
        delta = *sweep;
        *sweep = cur_depth;
        cur_depth += delta;
    }

    return DepthStatus::Success(depth_of_region);
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
