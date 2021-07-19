use crate::MEGA_BYTE_SIZE;
use crate::{Codecs, Writer};
use bam_tools::record::bamrawrecord::BAMRawRecord;
use bam_tools::sorting::sort;
use bam_tools::Reader;
use std::borrow::Cow;
use std::fs::File;
use std::io::{BufReader, BufWriter};

const MEM_LIMIT: usize = 2000 * MEGA_BYTE_SIZE;

/// Converts BAM file to GBAM file. This uses the `bam_parallel`
/// reader.
pub fn bam_to_gbam(in_path: &str, out_path: &str, codec: Codecs) {
    let (mut bgzf_reader, mut writer) = get_bam_reader_gbam_writer(in_path, out_path, codec);

    bgzf_reader.read_header().unwrap();
    bgzf_reader.consume_reference_sequences().unwrap();

    let mut records = bgzf_reader.records();
    while let Some(Ok(rec)) = records.next_rec() {
        let wrapper = BAMRawRecord(Cow::Borrowed(rec));
        writer.push_record(&wrapper);
    }

    writer.finish().unwrap();
}

pub fn bam_sort_to_gbam(in_path: &str, out_path: &str, codec: Codecs) {
    let fin = File::open(in_path).expect("failed");
    let fout = File::create(out_path).expect("failed");

    let buf_reader = BufReader::new(fin);
    let buf_writer = BufWriter::new(fout);

    let mut writer = Writer::new(buf_writer, codec, 8);

    let tmp_dir_path = std::env::temp_dir();

    sort::sort_bam(
        MEM_LIMIT,
        buf_reader,
        &mut writer,
        tmp_dir_path,
        0,
        20,
        20,
        false,
        sort::SortBy::NameAndMatchMates,
    )
    .unwrap();

    writer.finish().unwrap();
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

    let bgzf_reader = Reader::new(buf_reader, 20);
    let writer = Writer::new(buf_writer, codec, 8);

    (bgzf_reader, writer)
}
