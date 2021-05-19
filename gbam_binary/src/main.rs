// use gbam_tools::bam_to_gbam;
use gbam_tools::{bam_to_gbam, Codecs, Fields, ParsingTemplate, Reader};
use std::path::PathBuf;
use std::time::Instant;
use structopt::StructOpt;

#[derive(StructOpt)]
struct Cli {
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

fn main() {
    let args = Cli::from_args();
    if args.convert {
        bam_to_gbam(
            args.in_path.as_path().to_str().unwrap(),
            args.out_path
                .expect("Output path is mandatory for this operation.")
                .as_path()
                .to_str()
                .unwrap(),
            Codecs::Lz4,
        );
    } else if args.test {
        let mut tmplt = ParsingTemplate::new();
        tmplt.set(&Fields::Pos, true);
        tmplt.set(&Fields::RawCigar, true);
        let mut reader = Reader::new_for_file(args.in_path.as_path().to_str().unwrap(), tmplt);
        let now = Instant::now();
        #[allow(unused_variables)]
        while let Some(rec) = reader.next_rec() {}
        println!(
            "GBAM. Time elapsed querying POS and RAWCIGAR field throughout whole file: {}ms",
            now.elapsed().as_millis()
        );
    }
}

// fn main() {
//     let args: Vec<String> = std::env::args().collect();
//     if (args.len() != 4) {
//         eprintln!("Wrong number of arguments.");
//         eprintln!("GBAM writer.\nUsage: ./gbam 'path to bam file' 'path to output gbam directory' 'Compression type [ZSTD, FLATE2]'\n");
//         return;
//     }
//     let bam_path = &args[1];
//     let dir_path = &args[2];
//     let compr_type = match args[3].as_str() {
//         "ZSTD" => Compression::ZSTD,
//         "FLATE2" => Compression::FLATE2,
//         _ => {
//             eprintln!("Unknown option.");
//             return;
//         }
//     };

//     let mut writers = Vec::<std::io::BufWriter<File>>::new();
//     for field in Fields::iterator() {
//         let cur_name: &str = &(*field as usize).to_string();
//         writers.push(std::io::BufWriter::new(
//             File::create(dir_path.to_owned() + cur_name).expect("Unable to create file"),
//         ));
//     }

//     let mut gbam_writer = Writer::new(writers, compr_type).unwrap();
//     test_noodles_reader(bam_path, &mut gbam_writer);
// }

// fn test_noodles_reader(path: &str, writer: &mut Writer) {
//     let f = File::open(path).expect("failed");
//     let buf_reader = std::io::BufReader::new(f);
//     let mut reader = ParallelReader::new(buf_reader, 100);
//     // let mut reader = Reader::new(File::open(path).expect("failed"));
//     // let mut reader = Reader::new(buf_reader);

//     let mut buf = RawRecord::from(Vec::<u8>::new());
//     let mut cur_val: usize = 0;

//     let mut magic = [0; 4];
//     reader.read_exact(&mut magic).expect("Failed to read.");

//     if &magic[..] != b"BAM\x01" {
//         panic!("invalid BAM header");
//     }

//     let l_text = reader.read_u32::<LittleEndian>().expect("Failed to read.");
//     let mut c_text = vec![0; l_text as usize];
//     reader.read_exact(&mut c_text).expect("Failed to read.");

//     let n_ref = reader.read_u32::<LittleEndian>().expect("Failed to read.");

//     for _ in 0..n_ref {
//         let l_name = reader.read_u32::<LittleEndian>().expect("Failed to read.");
//         // eprintln!("{}", n_ref);
//         let mut c_name = vec![0; l_name as usize];
//         reader.read_exact(&mut c_name).expect("Failed to read.");
//         reader.read_u32::<LittleEndian>().expect("Failed to read.");
//     }

//     loop {
//         let block_size = match reader.read_u32::<LittleEndian>() {
//             Ok(bs) => bs as usize,
//             Err(ref e) if e.kind() == std::io::ErrorKind::UnexpectedEof => 0,
//             Err(e) => panic!(e),
//         };
//         if block_size == 0 {
//             eprintln!("{}", cur_val);
//             // EOF
//             return;
//         }
//         cur_val += block_size;
//         // println!("STATUS: {}\n", cur_val);

//         buf.resize(block_size);
//         reader.read_exact(&mut buf).expect("FAILED TO READ.");
//         // write!(writer, "{:?}", buf).expect("Output failed");
//         writer.write_record(&buf);
//     }
// }

//

// Store all fixed fields before variable sized ones.
// Store jump table before variable sized fields.

//
