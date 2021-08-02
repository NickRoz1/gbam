use std::collections::HashSet;

use bam_tools::record::fields::Fields;

use crate::reader::{reader::Reader, record::GbamRecord};

/// This module provides function for fast querying of read depth.

/// Get depth at position. The file should be sorted!
pub fn get_depth(reader: &mut Reader, chr_num: usize, pos: usize) -> Option<usize> {
    if !reader
        .parsing_template
        .check_if_active(&[Fields::RefID, Fields::Pos, Fields::RawSeqLen])
    {
        panic!("The reader should have parsing template which includes REFID, POS and L_SEQ.");
    }
}

fn find_refid(reader: &mut Reader, chr_num: usize, buf: &mut GbamRecord) -> Option<usize> {
    let mut l = 0;
    let mut r = reader.rec_num;

    while l <= r {
        let mid = (l + r) / 2;
        reader.fill_record(mid, buf);
        if buf.refid < chr_num {
            l = mid + 1;
        } else {
            r = mid - 1;
        }
    }

    match buf.refid == chr_num {
        true => Some(l),
        false => None,
    }
}
