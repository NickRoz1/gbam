// use gbam_tools::bam_to_gbam;
use bam_tools::record::fields::Fields;
use gbam_tools::{
    bam::bam_to_gbam::bam_sort_to_gbam,
    bam::gbam_to_bam::gbam_to_bam,
    query::depth::main_depth,
    reader::{parse_tmplt::ParsingTemplate, reader::Reader},
    {bam_to_gbam, Codecs},
    query::flagstat::collect_stats,
};


use std::path::PathBuf;
use std::time::Instant;
use std::fs::File;
use structopt::StructOpt;

#[derive(StructOpt)]
struct Cli {
    /// Sort BAM file before converting it to GBAM.
    #[structopt(short, long)]
    sort: bool,
    /// Specify which kind of temporary medium to use while sorting: ram, lz4_ram, file, lz4_file
    #[structopt(long)]
    sort_temp_mode: Option<String>,
    /// Determines whether conversion is requested
    #[structopt(short, long)]
    convert_to_gbam: bool,
    /// Convert to bam
    #[structopt(long)]
    convert_to_bam: bool,
    /// Perform the test
    #[structopt(short, long)]
    test: bool,
    /// Get depth at position.
    #[structopt(short, long)]
    depth: bool,
    /// Collect statistic from flag field from all records in the file.
    #[structopt(short, long)]
    flagstat: bool,
    /// The path to the BAM file to read
    #[structopt(parse(from_os_str))]
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
    /// Depth query. Filter reads with map quality lower than.
    #[structopt(long)]
    mapq: Option<u32>,
    /// Depth query. Number of threads to use. WARNING: each thread will attempt to allocate up to 1GB.
    #[structopt(long)]
    thread_num: Option<usize>,
    /// Sort temp directory.
    #[structopt(long, parse(from_os_str))]
    temp_dir: Option<PathBuf>,
}

/// Limited wrapper of `gbam_tools` converts BAM file to GBAM
/// file. Also limited tests may be run.
fn main() {
    let args = Cli::from_args();
    if args.convert_to_gbam {
        convert(args);
    } else if args.test {
        test(args);
    } else if args.depth {
        depth(args);
    } else if args.convert_to_bam {
        convert_to_bam(args);
    } else if args.flagstat {
        flagstat(args);
    }
}

fn convert(args: Cli) {
    let in_path = args
        .in_path
        .as_path()
        .to_str()
        .expect("Couldn't parse input path");
    let out_path = args
        .out_path
        .as_ref()
        .expect("Output path is mandatory for this operation.")
        .as_path()
        .to_str()
        .unwrap();
    if args.sort {
        bam_sort_to_gbam(in_path, out_path, Codecs::Lz4, args.sort_temp_mode, args.temp_dir);
    } else {
        bam_to_gbam(in_path, out_path, Codecs::Lz4);
    }
}

fn convert_to_bam(args: Cli) {
    let in_path = args
        .in_path
        .as_path()
        .to_str()
        .expect("Couldn't parse input path.");
    let out_path = args
        .out_path
        .as_ref()
        .expect("Output path is mandatory for this operation.")
        .as_path()
        .to_str()
        .unwrap();
    gbam_to_bam(in_path, out_path);
}

fn flagstat(args: Cli) {
    let in_path = args
        .in_path
        .as_path()
        .to_str()
        .expect("Couldn't parse input path.");

    let mut tmplt = ParsingTemplate::new();
    tmplt.set(&Fields::Flags, true);
    tmplt.set(&Fields::RefID, true);
    tmplt.set(&Fields::NextRefID, true);
    tmplt.set(&Fields::Mapq, true);

    let file = File::open(in_path).unwrap();

    collect_stats(file);
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
    #[allow(unused_mut)]
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
    drop(records);
}

fn depth(args: Cli) {
    let in_path = args.in_path.as_path().to_str().unwrap();
    let gbam_file = File::open(in_path).unwrap();
    main_depth(gbam_file, args.bed_file.as_ref(), args.query, args.mapq, args.out_path, args.thread_num);
}