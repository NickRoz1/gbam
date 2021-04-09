use crate::{Fields, RawRecord, Writer};
use byteorder::{ByteOrder, LittleEndian, ReadBytesExt, WriteBytesExt};
use noodles_bgzf::{ParallelReader, Reader};
use std::fs::File;
use std::io::{Read, Seek, Write};

use libc;

use libc::c_char;
use std::ffi::CStr;
use std::str;

/// Wrapper for Python FFI
#[no_mangle]
pub extern "C" fn bam_to_gbam_python(str1: *const c_char, str2: *const c_char) {
    let buf1 = unsafe { CStr::from_ptr(str1).to_bytes() };
    let buf2 = unsafe { CStr::from_ptr(str2).to_bytes() };
    let string1 = str::from_utf8(buf1).unwrap().to_owned();
    let string2 = str::from_utf8(buf2).unwrap().to_owned();
    bam_to_gbam(string1, string2);
}

/// Converts BAM file to GBAM file
pub fn bam_to_gbam(in_path: String, out_path: String) {
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
    let mut writer = Writer::new(buf_writer);

    let mut buf = RawRecord::from(Vec::<u8>::new());
    let mut cur_val: usize = 0;

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
            Err(e) => panic!(e),
        };
        if block_size == 0 {
            eprintln!("PARSED IN TOTAL: {}", cur_val);
            // EOF
            writer.finish();
            return ();
        }
        cur_val += block_size;
        // println!("STATUS: {}\n", cur_val);

        buf.resize(block_size);
        reader.read_exact(&mut buf).expect("FAILED TO READ.");
        // write!(writer, "{:?}", buf).expect("Output failed");
        writer.push_record(&buf);
    }
}
