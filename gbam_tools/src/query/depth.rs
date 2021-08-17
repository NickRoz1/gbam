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
    find_refid(reader, chr, buf)
}

//         https://github.com/biod/sambamba/blob/3eff9a2d8bb3097b92c72752be3c6b42dd1c59b7/BioD/bio/std/hts/bam/read.d#L902

fn find_refid(reader: &mut Reader, chr: String, buf: &mut GbamRecord) -> Option<usize> {
    let mut l = 0;
    let mut r = reader.rec_num;

    let chr_num = match reader.file_meta.get_ref_id(chr) {
        Some(num) => num,
        None => return None,
    };
    while l <= r {
        let mid = (l + r) / 2;
        reader.fill_record(mid, buf);
        if buf.refid.unwrap() < chr_num {
            l = mid + 1;
        } else {
            r = mid - 1;
        }
    }

    match buf.refid.unwrap() == chr_num {
        true => Some(l),
        false => None,
    }
}

/// Implement custom is_consuming_reference cigar parser. No need to do
/// LittleEndian read, when the only thing we need is 1 bytes (4 bits to be
/// precise).
/// Then, I guess the stuff from sambamba should do the trick for the rest.
pub fn l() {}
