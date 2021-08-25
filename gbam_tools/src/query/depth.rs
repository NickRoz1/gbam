use std::{cmp::Ordering, collections::HashSet, convert::TryFrom, time::Instant};

use bam_tools::record::fields::Fields;

use crate::reader::{reader::generate_block_treemap, reader::Reader, record::GbamRecord};

use crate::reader::record::copying;

/// This module provides function for fast querying of read depth.

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

    // Used to quickly determine what block record belongs to.
    let blocks = generate_block_treemap(&reader.file_meta, &Fields::RawCigar);

    let now = Instant::now();

    // Find number of record when records with this ref_id begin.
    if let Some(region_start) = find_refid(reader, chr_num) {
        // Find number of record when records with next ref_id start, or the
        // total amount of records if there is no such chr.
        let region_end = match find_refid(reader, chr_num + 1) {
            Some(end) => end,
            None => reader.amount,
        };
        let res = Some(calc_depth_at_pos(reader, region_start, region_end, pos));
        unsafe {
            println!(
                "Spent: {} ms\nOf which spent copying: {}\n",
                now.elapsed().as_millis(),
                copying.as_millis()
            );
        }

        res
    } else {
        None
    }
}

fn find_refid(reader: &mut Reader, chr_num: i32) -> Option<usize> {
    let rec_num = reader.amount;
    let pred = |num: usize, buf: &mut GbamRecord| {
        reader.fill_record(num, buf);
        buf.refid.unwrap().cmp(&chr_num)
    };

    binary_search(0, rec_num, pred)
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

pub fn calc_depth_at_pos(
    reader: &mut Reader,
    mut region_start: usize,
    mut region_end: usize,
    pos: i32,
) -> usize {
    let mut buf = GbamRecord::default();
    let mut depth = 0;

    loop {
        if region_start >= reader.amount {
            return depth;
        }
        reader.fill_record(region_start, &mut buf);
        region_start += 1;

        let read_pos = buf.pos.unwrap();
        let base_cov = buf.cigar.as_ref().unwrap().base_coverage();

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

        if read_pos >= pos {
            return depth;
        }

        // https://github.com/biod/sambamba/blob/c795656721b3608ffe7765b6ab98502426d14131/BioD/bio/std/hts/bam/randomaccessmanager.d#L448
        if read_pos + i32::try_from(base_cov).unwrap() > pos {
            depth += 1;
        }
    }
}
