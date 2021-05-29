use crate::{Codecs, BAMRawRecord, Writer};
use bam_parallel::ParallelReader;
use byteorder::{LittleEndian, ReadBytesExt};
use std::fs::File;
use std::io::Read;

/// Converts BAM file to GBAM file. This uses the `bam_parallel`
/// reader.
pub fn bam_to_gbam(in_path: &str, out_path: &str, codec: Codecs) {
    // println!(
    //     "PYTHON FFI IS WORKING. INPUT PARAMETERS1 ARE: {} | {}",
    //     in_path, out_path
    // );
    // return ();
    let fin = File::open(in_path).expect("failed");
    let fout = File::create(out_path).expect("failed");

    let buf_reader = std::io::BufReader::new(fin);
    let buf_writer = std::io::BufWriter::new(fout);

    let mut reader = ParallelReader::new(buf_reader, 10);
    let mut writer = Writer::new(buf_writer, codec, 8);

    // Check for BAM magic header
    let mut buf = BAMRawRecord::from(Vec::<u8>::new());

    let mut magic = [0; 4];
    reader.read_exact(&mut magic).expect("Failed to read.");

    if &magic[..] != b"BAM\x01" {
        panic!("invalid BAM header");
    }

    let l_text = reader.read_u32::<LittleEndian>().expect("Failed to read.");
    let mut c_text = vec![0; l_text as usize];
    reader.read_exact(&mut c_text).expect("Failed to read.");

    let n_ref = reader.read_u32::<LittleEndian>().expect("Failed to read.");
    // eprintln!("{}", n_ref);

    for _ in 0..n_ref {
        let l_name = reader.read_u32::<LittleEndian>().expect("Failed to read.");
        let mut c_name = vec![0; l_name as usize];
        reader.read_exact(&mut c_name).expect("Failed to read.");
        reader.read_u32::<LittleEndian>().expect("Failed to read.");
    }

    loop {
        let block_size = match reader.read_u32::<LittleEndian>() {
            Ok(bs) => bs as usize,
            Err(ref e) if e.kind() == std::io::ErrorKind::UnexpectedEof => 0,
            Err(e) => panic!("{}", e),
        };
        if block_size == 0 {
            // EOF
            writer.finish().unwrap();
            return;
        }

        // Fetch the BAM read and push to the GBAM writer
        buf.resize(block_size);
        reader.read_exact(&mut buf).expect("FAILED TO READ.");
        // write!(writer, "{:?}", buf).expect("Output failed");
        writer.push_record(&buf);
    }
}
