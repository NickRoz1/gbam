use crate::MEGA_BYTE_SIZE;
use crate::{Codecs, Writer};
use bam_tools::parse_reference_sequences;
use bam_tools::record::bamrawrecord::BAMRawRecord;
use bam_tools::record::fields::{Fields, FIELDS_NUM};
use bam_tools::sorting::sort;
use bam_tools::Reader;
use std::borrow::Cow;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufWriter};

const MEM_LIMIT: usize = 2000 * MEGA_BYTE_SIZE;

/// Converts BAM file to GBAM file. This uses the `bam_parallel` reader.
pub fn bam_to_gbam(in_path: &str, out_path: &str, codec: Codecs) {
    let (mut bam_reader, mut writer) = get_bam_reader_gbam_writer(in_path, out_path, codec);

    let mut records = bam_reader.records();
    while let Some(Ok(rec)) = records.next_rec() {
        let wrapper = BAMRawRecord(Cow::Borrowed(rec));
        writer.push_record(&wrapper);
    }

    writer.finish().unwrap();
}

/// Converts BAM file to GBAM file. Sorts BAM file in process. This uses the `bam_parallel` reader.
pub fn bam_sort_to_gbam(in_path: &str, out_path: &str, codec: Codecs) {
    let fin_for_ref_seqs = File::open(in_path).expect("failed");
    let mut reader_for_header_only = Reader::new(fin_for_ref_seqs, 1);
    let (sam_header, ref_seqs, ref_seq_offset) =
        read_sam_header_and_ref_seqs(&mut reader_for_header_only);
    let only_text = &sam_header[..ref_seq_offset];

    let fin = File::open(in_path).expect("failed");
    let fout = File::create(out_path).expect("failed");

    let buf_reader = BufReader::new(fin);
    let buf_writer = BufWriter::new(fout);

    let mut writer = Writer::new(
        buf_writer,
        vec![codec; FIELDS_NUM],
        8,
        vec![Fields::RefID],
        ref_seqs,
        Vec::from(only_text),
    );

    let tmp_dir_path = std::env::temp_dir();

    sort::sort_bam(
        MEM_LIMIT,
        buf_reader,
        &mut writer,
        tmp_dir_path,
        0,
        1,
        20,
        false,
        sort::SortBy::NameAndMatchMates,
    )
    .unwrap();

    writer.finish().unwrap();
}

/// Consumes SAM header from input BAM reader.
///
///
/// # Returns
/// Tuple of 3 elements.
///
/// **tuple.0** -> all bytes of BAM header including reference sequences.
///
/// **tuple.1** -> parsed reference sequences from BAM header.
///
/// **tuple.2** -> offset to reference sequences in tuple.0. It's before n_ref uint32_t.
fn read_sam_header_and_ref_seqs(reader: &mut Reader) -> (Vec<u8>, Vec<(String, u32)>, usize) {
    let (bytes_of_header, ref_sequences_offset) = reader.read_header().unwrap();
    let sequences = parse_reference_sequences(&bytes_of_header[ref_sequences_offset..]).unwrap();
    (bytes_of_header, sequences, ref_sequences_offset)
}

fn get_bam_reader_gbam_writer(
    in_path: &str,
    out_path: &str,
    codec: Codecs,
) -> (Reader, Writer<BufWriter<File>>) {
    let fin = File::open(in_path).expect("failed");
    let fout = File::create(out_path).expect("failed");

    let buf_reader = BufReader::new(fin);
    let buf_writer = BufWriter::new(fout);

    let mut bgzf_reader = Reader::new(buf_reader, 1);

    let (sam_header, ref_seqs, ref_seq_offset) = read_sam_header_and_ref_seqs(&mut bgzf_reader);

    let only_text = &sam_header[..ref_seq_offset];

    let writer = Writer::new(
        buf_writer,
        vec![codec; FIELDS_NUM],
        8,
        vec![Fields::RefID],
        ref_seqs,
        Vec::from(only_text),
    );

    (bgzf_reader, writer)
}
