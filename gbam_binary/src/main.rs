// use gbam_tools::bam_to_gbam;
use bam_tools::record::fields::Fields;
use gbam_tools::{
    bam::bam_to_gbam::bam_sort_to_gbam,
    query::depth::{get_region_depth, DepthStatus},
    reader::{parse_tmplt::ParsingTemplate, reader::Reader, records::Records},
    utils::bed,
    {bam_to_gbam, Codecs},
};

use std::fs::File;
use std::path::PathBuf;
use std::time::Instant;
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
    #[structopt(parse(from_os_str))]
    in_path: PathBuf,
    /// The path to write output GBAM file
    #[structopt(parse(from_os_str))]
    out_path: Option<PathBuf>,
    /// Depth query. Example: chr1:54-54, or chrX:1258-9999
    #[structopt(short, long)]
    query: Option<String>,
    /// Depth query. Example: chr1:54, or chrX:1258
    #[structopt(parse(from_os_str))]
    bed_file: Option<PathBuf>,
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
    // let now = Instant::now();
    // println!(
    //     "GBAM. Time elapsed querying depth {}ms",
    //     now.elapsed().as_millis()
    // );
    let mut buf = Vec::new();
    if !queries.is_empty() {
        println!("REFID\tPOS\tDEPTH");
    }
    for query in queries.iter() {
        match get_region_depth(&mut reader, query, buf) {
            DepthStatus::Success(region_depth) => {
                for (i, depth) in region_depth.iter().enumerate() {
                    // Start pos + current offset.
                    let position = query.1 as usize + i;
                    println!("{}\t{}\t{}", query.0, position, depth);
                }
                buf = region_depth;
            }
            DepthStatus::None(buffer) => {
                buf = buffer;
                println!(
                    "Unable to locate ref id {} or range {}-{}",
                    query.0, query.1, query.2
                );
            }
        }
    }
}

fn panic_err() {
    panic!("The query you entered is incorrect. The format is as following: <ref name>:<position>\ne.g. chr1:1257\n");
}
