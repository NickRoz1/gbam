// use gbam_tools::bam_to_gbam;
use bam_tools::record::fields::Fields;
use gbam_tools::bam::bam_to_gbam::bam_sort_to_gbam;
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
    /// The path to the BAM file to read
    #[structopt(parse(from_os_str))]
    in_path: PathBuf,
    /// The path to write output GBAM file
    #[structopt(parse(from_os_str))]
    out_path: Option<PathBuf>,
}

/// Limited wrapper of `gbam_tools` converts BAM file to GBAM
/// file. Also limited tests may be run.
fn main() {
    let args = Cli::from_args();
    if args.convert {
        convert(args);
    } else if args.test {
        test(args);
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
    tmplt.set(&Fields::Pos, true);
    tmplt.set(&Fields::RawCigar, true);
    let file = File::open(args.in_path.as_path().to_str().unwrap()).unwrap();
    let mut reader = Reader::new(file, tmplt).unwrap();
    let mut records = reader.records();
    let now = Instant::now();
    let mut u = 1;
    #[allow(unused_variables)]
    while let Some(rec) = records.next_rec() {
        u += 1;
    }
    println!(
        "GBAM. Time elapsed querying POS and RAWCIGAR field throughout whole file: {}ms",
        now.elapsed().as_millis()
    );
}
