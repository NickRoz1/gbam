use crate::stats::{refid_comparator, StatsComparator};
use crate::MEGA_BYTE_SIZE;
use crate::{Codecs, Writer};
use bam_tools::record::bamrawrecord::BAMRawRecord;
use bam_tools::record::fields::{Fields, FIELDS_NUM};
use bam_tools::sorting::sort;
use rust_htslib::bam::{FetchDefinition, IndexedReader, Read, Record};
use rust_htslib::errors::Error;
use std::borrow::Cow;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufWriter};

const MEM_LIMIT: usize = 2000 * MEGA_BYTE_SIZE;

/// Converts BAM file to GBAM file. This uses the `bam_parallel`
/// reader.
pub fn bam_to_gbam(in_path: &str, out_path: &str, codec: Codecs) {
    //  let mut bam = IndexedReader::from_path(in_path).unwrap();

    //  bam.fetch(FetchDefinition::All);

    //  let mut record = Record::new();
    //  while let Some(result) = bam.read(&mut record) {
    //      match result {
    //          Ok(_) => {
    //              println!("Read sequence: {:?}", record.seq().as_bytes());
    //          },
    //          Err(_) => panic!("BAM parsing failed...")
    //      }
    //  }

    // let fout = File::create(out_path).expect("failed");
    // let buf_writer = BufWriter::new(fout);

    // let mut stats_collectors = HashMap::<Fields, StatsComparator>::new();
    // stats_collectors.insert(Fields::RefID, refid_comparator);

    // let writer = Writer::new(
    //     buf_writer,
    //     vec![codec; FIELDS_NUM],
    //     8,
    //     stats_collectors,

    // );

    // writer

    // let buf_reader = BufReader::new(
    // let mut records = bgzf_reader.records();
    // while let Some(Ok(rec)) = records.next_rec() {
    //     let wrapper = BAMRawRecord(Cow::Borrowed(rec));
    //     writer.push_record(&wrapper);
    // }

    // writer.finish().unwrap();
}

pub fn bam_sort_to_gbam(in_path: &str, out_path: &str, codec: Codecs) {
    // let fin_for_ref_seqs = File::open(in_path).expect("failed");
    // let ref_seqs = extract_ref_seqs(BufReader::new(fin_for_ref_seqs));

    // let fin = File::open(in_path).expect("failed");
    // let fout = File::create(out_path).expect("failed");

    // let buf_reader = BufReader::new(fin);
    // let buf_writer = BufWriter::new(fout);

    // let mut stats_collectors = HashMap::<Fields, StatsComparator>::new();
    // stats_collectors.insert(Fields::RefID, refid_comparator);

    // let mut writer = Writer::new(
    //     buf_writer,
    //     vec![codec; FIELDS_NUM],
    //     8,
    //     stats_collectors,
    //     ref_seqs,
    // );

    // let tmp_dir_path = std::env::temp_dir();

    // sort::sort_bam(
    //     MEM_LIMIT,
    //     buf_reader,
    //     &mut writer,
    //     tmp_dir_path,
    //     0,
    //     20,
    //     20,
    //     false,
    //     sort::SortBy::NameAndMatchMates,
    // )
    // .unwrap();

    // writer.finish().unwrap();
}

// fn extract_ref_seqs(reader: BufReader<File>) -> Vec<(String, i32)> {
//     // let mut parallel_reader = Reader::new(reader, 1);
//     // parallel_reader.read_header().unwrap();
//     // parallel_reader.parse_reference_sequences().unwrap()
// }

// fn get_bam_reader_gbam_writer(out_path: &str, codec: Codecs) -> Writer<BufWriter<File>> {

// }
