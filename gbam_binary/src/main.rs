// use gbam_tools::bam_to_gbam;
use bam_tools::record::fields::Fields;
use gbam_tools::{
    bam::bam_to_gbam::bam_sort_to_gbam,
    query::depth::get_regions_depths,
    reader::{parse_tmplt::ParsingTemplate, reader::Reader, records::Records},
    utils::bed,
    {bam_to_gbam, Codecs},
};

use bam_tools::record::bamrawrecord::BAMRawRecord;
use rust_htslib::bam::{FetchDefinition, Read, Record};
use rust_htslib::errors::Error;
use std::io::{BufReader, BufWriter};
use std::path::PathBuf;
use std::time::Instant;
use std::{borrow::Cow, fs::File};
use structopt::StructOpt;

#[derive(StructOpt)]
struct Cli {
    /// Sort BAM file before converting it to GBAM.
    #[structopt(short, long)]
    sort: bool,
    /// Determines whether conversion is requested
    #[structopt(short, long)]
    convert: bool,
    /// Perform the test
    #[structopt(short, long)]
    test: bool,
    /// Get depth at position.
    #[structopt(short, long)]
    depth: bool,
    /// The path to the BAM file to read
    #[structopt(short, parse(from_os_str))]
    in_path: PathBuf,
    /// The path to write output GBAM file
    #[structopt(short, parse(from_os_str))]
    out_path: Option<PathBuf>,
    /// Depth query. Example: chr1:54-54, or chrX:1258-9999
    #[structopt(short, long)]
    query: Option<String>,
    /// Depth query. Example: chr1:54, or chrX:1258
    #[structopt(short, parse(from_os_str))]
    bed_file: Option<PathBuf>,
    /// Test speed of reading whole BAM file with my reader.
    #[structopt(long)]
    test_speed_my_reader: bool,
    /// Test speed of reading whole BAM file with rust hts reader.
    #[structopt(long)]
    test_speed_rust_hts: bool,
}

/// Limited wrapper of `gbam_tools` converts BAM file to GBAM
/// file. Also limited tests may be run.
fn main() {
    let args = Cli::from_args();
    if args.convert {
        convert(args);
    } else if args.test {
        test(args);
    } else if args.depth {
        depth(args);
    } else if args.test_speed_my_reader {
        fn_test_read_time(Kind::my_reader, args.in_path.as_path().to_str().unwrap());
    } else if args.test_speed_rust_hts {
        fn_test_read_time(Kind::rust_hts, args.in_path.as_path().to_str().unwrap());
    }
}

enum Kind {
    my_reader,
    rust_hts,
}

fn fn_test_read_time(kind: Kind, in_path: &str) {
    match kind {
        Kind::my_reader => {
            let fin = File::open(in_path).expect("failed");

            let now = Instant::now();
            let buf_reader = BufReader::new(fin);

            let mut bgzf_reader = bam_tools::Reader::new(buf_reader, 20);

            bgzf_reader.read_header().unwrap();
            bgzf_reader.parse_reference_sequences().unwrap();

            let mut records = bgzf_reader.records();
            let mut counter = 0;
            while let Some(Ok(rec)) = records.next_rec() {
                let wrapper = BAMRawRecord(Cow::Borrowed(rec));
                counter += 1;
            }

            let cur = now.elapsed().as_millis();
            println!("My reader took : {}ms, to read this file.", cur);
        }
        Kind::rust_hts => {
            let now = Instant::now();

            let mut bam = rust_htslib::bam::Reader::from_path(in_path).unwrap();

            let pool = rust_htslib::tpool::ThreadPool::new(20).unwrap();

            bam.set_thread_pool(&pool).unwrap();

            let mut record = Record::new();
            let mut counter = 0;
            while let Some(result) = bam.read(&mut record) {
                match result {
                    Ok(_) => {
                        counter += 1;
                        // println!("Read sequence: {:?}", record.seq().as_bytes());
                    }
                    Err(_) => panic!("BAM parsing failed..."),
                }
            }
            let cur = now.elapsed().as_millis();
            println!("Rust_hts took : {}ms, to read this file.", cur);
        }
    }
}

fn convert(args: Cli) {
    let in_path = args.in_path.as_path().to_str().unwrap();
    let out_path = args
        .out_path
        .as_ref()
        .expect("Output path is mandatory for this operation.")
        .as_path()
        .to_str()
        .unwrap();
    if args.sort {
        bam_sort_to_gbam(in_path, out_path, Codecs::Lz4)
    } else {
        bam_to_gbam(in_path, out_path, Codecs::Lz4);
    }
}

fn test(args: Cli) {
    let mut tmplt = ParsingTemplate::new();
    tmplt.set(&Fields::RefID, true);
    tmplt.set(&Fields::Pos, true);
    tmplt.set(&Fields::RawCigar, true);
    let file = File::open(args.in_path.as_path().to_str().unwrap()).unwrap();
    let mut reader = Reader::new(file, tmplt).unwrap();
    let mut records = reader.records();
    let now = Instant::now();
    let mut u = 0;
    #[allow(unused_variables)]
    while let Some(rec) = records.next_rec() {
        // u += rec.cigar.as_ref().unwrap().as_bytes().len();
    }
    println!("Record count {}", u);
    println!(
        "GBAM. Time elapsed querying POS and RAWCIGAR field throughout whole file: {}ms",
        now.elapsed().as_millis()
    );
}

fn depth(args: Cli) {
    let mut tmplt = ParsingTemplate::new();
    tmplt.set(&Fields::RefID, true);
    tmplt.set(&Fields::Pos, true);
    tmplt.set(&Fields::RawCigar, true);

    let in_path = args.in_path.as_path().to_str().unwrap();
    let file = File::open(in_path).unwrap();
    let mut reader = Reader::new(file, tmplt).unwrap();
    let mut queries = Vec::new();
    if let Some(bed_path) = args.bed_file {
        queries = bed::parse_bed_from_file(&bed_path).expect("BED file is corrupted.");
    } else {
        queries.push(
            bed::parse_region_query_owned(&args.query.unwrap())
                .ok()
                .ok_or_else(|| panic_err())
                .unwrap(),
        );
    }

    if !queries.is_empty() {
        println!("REFID\tPOS\tDEPTH");
    }
    let now = Instant::now();
    get_regions_depths(&mut reader, &queries);
    println!(
        "GBAM. Time elapsed querying depth {}ms",
        now.elapsed().as_millis()
    );
}

fn panic_err() {
    panic!("The query you entered is incorrect. The format is as following: <ref name>:<position>\ne.g. chr1:1257\n");
}
