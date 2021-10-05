// use gbam_tools::bam_to_gbam;
use bam_tools::record::fields::Fields;
use gbam_tools::bam::bam_to_gbam::bam_sort_to_gbam;
use gbam_tools::query::depth::get_depth;
use gbam_tools::reader::{parse_tmplt::ParsingTemplate, reader::Reader, records::Records};
use gbam_tools::{bam_to_gbam, Codecs};
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
    /// Depth query. Example: chr1:54, or chrX:1258
    #[structopt(short, long)]
    query: Option<String>,
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

    let (ref_id, pos) = parse_query(args.query.clone().unwrap());
    let now = Instant::now();
    if let Some(depth) = get_depth(&mut reader, ref_id, pos) {
        println!(
            "GBAM. Time elapsed querying depth {}ms",
            now.elapsed().as_millis()
        );
        println!(
            "The depth at <{}> is equal to: {}",
            args.query.unwrap(),
            depth
        );
        println!("test");
    } else {
        println!("No reads cover position {}", args.query.unwrap());
    }
}

fn parse_query(q: String) -> (String, i32) {
    let mut parts = q.split(':');
    if parts.clone().count() != 2 {
        panic_err();
    }
    let ref_id = parts.next().unwrap().to_owned();
    let pos = parts.next().unwrap().parse::<i32>().unwrap();

    (ref_id, pos)
}

fn panic_err() {
    panic!("The query you entered is incorrect. The format is as following: <ref name>:<position>\ne.g. chr1:1257\n");
}
