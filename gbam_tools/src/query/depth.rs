use std::{collections::HashSet, convert::TryFrom};

use bam_tools::record::fields::Fields;

use crate::reader::{reader::Reader, record::GbamRecord};

/// This module provides function for fast querying of read depth.

/// Get depth at position. The file should be sorted!
pub fn get_depth(reader: &mut Reader, chr: String, pos: i32) -> Option<usize> {
    if !reader
        .parsing_template
        .check_if_active(&[Fields::RefID, Fields::Pos, Fields::RawCigar])
    {
        panic!("The reader should have parsing template which includes REFID, POS and CIGAR.");
    }

    if let Some(region_start) = find_refid(reader, chr) {
        calc_depth_at_pos(reader, region_start, pos)
    } else {
        None
    }
}

fn find_refid(reader: &mut Reader, chr: String) -> Option<usize> {
    let mut l = 0;
    let mut r = reader.amount;

    let mut buf = GbamRecord::default();

    // https://github.com/biod/sambamba/blob/3eff9a2d8bb3097b92c72752be3c6b42dd1c59b7/BioD/bio/std/hts/bam/read.d#L902
    let chr_num = match reader.file_meta.get_ref_id(chr) {
        Some(num) => num,
        None => return None,
    };
    while l < r {
        let mid = (l + r) / 2;
        println!("Mid - {:?}", mid);
        reader.fill_record(mid, &mut buf);
        if buf.refid.unwrap() < chr_num {
            l = mid + 1;
        } else {
            r = mid;
        }
    }
    // reader.fill_record(0, &mut buf);
    println!("Hey {} - {}", chr_num, buf.refid.unwrap());
    match buf.refid.unwrap() == chr_num {
        true => Some(l),
        false => None,
    }
}

pub fn calc_depth_at_pos(reader: &mut Reader, mut region_start: usize, pos: i32) -> Option<usize> {
    let mut buf = GbamRecord::default();
    let mut depth = 0;

    loop {
        if region_start >= reader.amount {
            break;
        }
        reader.fill_record(region_start, &mut buf);
        region_start += 1;

        let read_pos = buf.pos.unwrap();
        let base_cov = buf.cigar.as_ref().unwrap().base_coverage();

        // println!(
        //     "REF_ID: {}\nPOS: {}\nbase_cov: {}\ncigar_len: {}\nsearch pos: {}\n",
        //     buf.refid.unwrap(),
        //     read_pos,
        //     base_cov,
        //     buf.cigar.as_ref().unwrap().0.len(),
        //     pos
        // );

        if read_pos > pos {
            return Some(depth);
        }

        // https://github.com/biod/sambamba/blob/c795656721b3608ffe7765b6ab98502426d14131/BioD/bio/std/hts/bam/randomaccessmanager.d#L448
        if read_pos + i32::try_from(base_cov).unwrap() > pos {
            depth += 1;
        }
    }

    None
}
