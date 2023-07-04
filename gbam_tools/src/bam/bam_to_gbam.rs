use crate::MEGA_BYTE_SIZE;
use crate::{Codecs, Writer};
use bam_tools::parse_reference_sequences;
use bam_tools::record::bamrawrecord::BAMRawRecord;
use bam_tools::record::fields::{Fields, FIELDS_NUM};
use bam_tools::sorting::sort;
use bam_tools::sorting::sort::TempFilesMode;
use bam_tools::Reader;
use std::borrow::Cow;
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::PathBuf;
use std::str::FromStr;
use tempdir::TempDir;


const MEM_LIMIT: usize = 2000 * MEGA_BYTE_SIZE;

/// Converts BAM file to GBAM file. This uses the `bam_parallel` reader.
pub fn bam_to_gbam(in_path: &str, out_path: &str, codec: Codecs, full_command: String) {
    let (mut bam_reader, mut writer) = get_bam_reader_gbam_writer(in_path, out_path, codec, full_command);

    let mut records = bam_reader.records();
    while let Some(Ok(rec)) = records.next_rec() {
        let wrapper = BAMRawRecord(Cow::Borrowed(rec));
        writer.push_record(&wrapper);
    }

    writer.finish().unwrap();
}

/// Converts BAM file to GBAM file. Sorts BAM file in process. This uses the `bam_parallel` reader.
pub fn bam_sort_to_gbam(in_path: &str, out_path: &str, codec: Codecs, mut sort_temp_mode: Option<String>, temp_dir: Option<PathBuf>, full_command: String) {
    let fin_for_ref_seqs = File::open(in_path).expect("failed");
    
    let mut reader_for_header_only = Reader::new(fin_for_ref_seqs, 1, None);
    let (sam_header, ref_seqs, ref_seq_offset) =
        read_sam_header_and_ref_seqs(&mut reader_for_header_only);
    let only_text = &sam_header[..ref_seq_offset];

    let fin = File::open(in_path).expect("failed");
    let fout = File::create(out_path).expect("failed");

    let file_size = fin.metadata().unwrap().len();

    let buf_reader = BufReader::new(fin);
    let buf_writer = BufWriter::new(fout);

    let mut writer = Writer::new(
        buf_writer,
        vec![codec; FIELDS_NUM],
        8,
        vec![Fields::RefID],
        ref_seqs,
        Vec::from(only_text),
        full_command,
        true
    );

    let tmp_dir_path = temp_dir.map_or(std::env::temp_dir(), |path| path);
    if sort_temp_mode.is_none() {
        sort_temp_mode = Some(String::from_str("file").unwrap());
    }
    let tmp_medium_mode = match sort_temp_mode.unwrap().as_str() {
        "file" => TempFilesMode::RegularFiles,
        "lz4_file" => TempFilesMode::LZ4CompressedFiles,
        "ram" => TempFilesMode::InMemoryBlocks,
        "lz4_ram" => TempFilesMode::InMemoryBlocksLZ4,
        _ => panic!("Unknown sort_temp_mode mode."),
    };

    let dir = TempDir::new_in(tmp_dir_path, "BAM sort temporary directory.").unwrap();

    sort::sort_bam(
        MEM_LIMIT,
        buf_reader,
        &mut writer,
        &dir,
        0,
        4,
        tmp_medium_mode,
        sort::SortBy::CoordinatesAndStrand,
        Some(file_size)
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
    full_command: String,
) -> (Reader, Writer<BufWriter<File>>) {
    let fin = File::open(in_path).expect("failed");
    let fout = File::create(out_path).expect("failed");

    let file_size = fin.metadata().unwrap().len();

    let buf_reader = BufReader::new(fin);
    let buf_writer = BufWriter::new(fout);

    let mut bgzf_reader = Reader::new(buf_reader, 4, Some(file_size));

    let (sam_header, ref_seqs, ref_seq_offset) = read_sam_header_and_ref_seqs(&mut bgzf_reader);

    let only_text = &sam_header[..ref_seq_offset];

    let writer = Writer::new(
        buf_writer,
        vec![codec; FIELDS_NUM],
        8,
        vec![Fields::RefID],
        ref_seqs,
        Vec::from(only_text),
        full_command,
        false,
    );

    (bgzf_reader, writer)
}
