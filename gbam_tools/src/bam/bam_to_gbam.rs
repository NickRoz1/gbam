use crate::MEGA_BYTE_SIZE;
use crate::{Codecs, Writer};
use crate::{SamRecord, SamHeader, SamRawRecords, SamReader};
use bam_tools::parse_reference_sequences;
use bam_tools::record::bamrawrecord::BAMRawRecord;
use bam_tools::record::fields::{Fields, FIELDS_NUM};
use bam_tools::sorting::sort;
use bam_tools::sorting::sort::TempFilesMode;
use bam_tools::Reader;
use std::borrow::Cow;
use std::fs::File;
use std::io::{self, Write, BufReader, BufWriter, Read, stdin};
use std::path::PathBuf;
use std::str::FromStr;
use tempdir::TempDir;

use byteorder::{LittleEndian, WriteBytesExt};
use std::convert::TryInto;


const MEM_LIMIT: usize = 2000 * MEGA_BYTE_SIZE;

/// Converts BAM file to GBAM file. This uses the `bam_parallel` reader.
pub fn bam_to_gbam(in_path: &str, out_path: &str, codec: Codecs, full_command: String) {
    let (mut bam_reader, mut writer) = get_bam_reader_gbam_writer(in_path, out_path, codec, full_command);

    let mut records = bam_reader.records();
    while let Some(Ok(rec)) = records.next_rec() {
        let wrapper: BAMRawRecord<'_> = BAMRawRecord(Cow::Borrowed(rec));
        writer.push_record(&wrapper);
    }

    writer.finish().unwrap();
}

pub fn sam_to_gbam(in_path: Option<&str>, out_path: &str, codec: Codecs, full_command: String) {
    let (mut sam_reader, mut writer) = get_sam_reader_gbam_writer(in_path, out_path, codec, full_command);
    let ref_names = sam_reader.reference_names();

    for result in sam_reader.records() {
        match result {
            Ok(sam_record) => {
                // Convert SAM to BAM record binary
                let bam_bytes = convert_sam_record_to_bam(&sam_record, &ref_names); // Your implementation
                // validate_bam_record(bam_bytes);
                let wrapper = BAMRawRecord::from(bam_bytes[4..].to_vec());
                writer.push_record(&wrapper);
            }
            Err(e) => {
                eprintln!("Skipping invalid record: {}", e);
                continue;
            }
        }
    }

    writer.finish().unwrap();
}

fn count_cigar_ops(cigar: &str) -> u16 {
    let mut count = 0;
    let mut num = false;

    for c in cigar.chars() {
        if c.is_ascii_digit() {
            num = true;
        } else if num {
            count += 1;
            num = false;
        } else {
            // CIGAR like "M" or "X" without number — invalid
            return 0;
        }
    }

    count
}

/// Converts a parsed SamRecord into a BAM binary record
pub fn convert_sam_record_to_bam(record: &SamRecord, ref_names: &[String]) -> Vec<u8> {
    let tid = ref_names
        .iter()
        .position(|r| r == &record.rname)
        .unwrap_or(-1_i32 as usize) as i32;

    let next_tid = match record.rnext.as_str() {
        "*" => -1,
        "=" => tid,
        name => ref_names
            .iter()
            .position(|r| r == name)
            .unwrap_or(-1_i32 as usize) as i32,
    };

    let l_read_name = record.qname.len() + 1; // including null terminator
    let n_cigar_op = count_cigar_ops(&record.cigar);
    let l_seq = record.seq.len();
    let bin = 0u16; // Optional: compute bin if needed

    let mut data = Vec::new();
    data.extend_from_slice(&[0u8; 4]); // placeholder for block_size

    // --- Fixed fields: total of 32 bytes ---
    data.write_i32::<LittleEndian>(tid).unwrap();                      // refID
    data.write_i32::<LittleEndian>(record.pos as i32 - 1).unwrap();   // pos (0-based)
    data.write_u8(l_read_name as u8).unwrap();                        // l_read_name
    data.write_u8(record.mapq).unwrap();                              // mapq
    data.write_u16::<LittleEndian>(bin).unwrap();                     // bin
    data.write_u16::<LittleEndian>(n_cigar_op).unwrap();              // n_cigar
    data.write_u16::<LittleEndian>(record.flag).unwrap();             // flag
    data.write_u32::<LittleEndian>(l_seq as u32).unwrap();            // l_seq
    data.write_i32::<LittleEndian>(next_tid).unwrap();                // next_refID
    data.write_i32::<LittleEndian>(record.pnext as i32 - 1).unwrap(); // next_pos
    data.write_i32::<LittleEndian>(record.tlen).unwrap();             // tlen

    // --- Variable-length fields ---
    data.extend_from_slice(record.qname.as_bytes()); // read_name
    data.push(0); // null terminator

    data.extend_from_slice(&encode_cigar(&record.cigar)); // CIGAR
    data.extend_from_slice(&encode_seq(&record.seq));     // SEQ
    data.extend_from_slice(&encode_qual(&record.qual));   // QUAL

    // Optional: Tags
    for tag_field in &record.tags {
        let parts: Vec<&str> = tag_field.splitn(3, ':').collect();
        if parts.len() != 3 {
            continue;
        }

        let tag = parts[0];
        let type_char = parts[1];
        let value = parts[2];

        match type_char {
            "A" => {
                data.extend_from_slice(tag.as_bytes());
                data.extend_from_slice(b"A");
                data.push(value.as_bytes()[0]);
            }
            "i" => {
                if let Ok(val) = value.parse::<i32>() {
                    data.extend_from_slice(tag.as_bytes());
                    data.extend_from_slice(b"i");
                    data.extend_from_slice(&val.to_le_bytes());
                }
            }
            "Z" => {
                data.extend_from_slice(tag.as_bytes());
                data.extend_from_slice(b"Z");
                data.extend_from_slice(value.as_bytes());
                data.push(0);
            }
            _ => {} // Skip unsupported types
        }
    }

    // Patch the block size at the beginning
    let block_size = (data.len() - 4) as u32;
    data[0..4].copy_from_slice(&block_size.to_le_bytes());

    data
}


fn encode_cigar(cigar: &str) -> Vec<u8> {
    let mut bytes = Vec::new();
    let mut num = String::new();

    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            num.push(ch);
        } else {
            let len: u32 = num.parse().unwrap_or(0);
            let op = match ch {
                'M' => 0,
                'I' => 1,
                'D' => 2,
                'N' => 3,
                'S' => 4,
                'H' => 5,
                'P' => 6,
                '=' => 7,
                'X' => 8,
                _ => panic!("Unsupported CIGAR op: {}", ch),
            };
            let code = (len << 4) | op;
            bytes.write_u32::<LittleEndian>(code).unwrap();
            num.clear();
        }
    }

    bytes
}

fn encode_seq(seq: &str) -> Vec<u8> {
    let mut result = Vec::new();
    let mut byte = 0u8;
    for (i, base) in seq.chars().enumerate() {
        let code = encode_base(base);
        if i % 2 == 0 {
            byte = code << 4;
        } else {
            byte |= code;
            result.push(byte);
            byte = 0;
        }
    }
    if seq.len() % 2 != 0 {
        result.push(byte);
    }
    result
}

fn encode_qual(qual: &str) -> Vec<u8> {
    if qual == "*" {
        vec![0xFF]
    } else {
        qual.bytes().map(|b| b.saturating_sub(33)).collect()
    }
}

fn encode_base(base: char) -> u8 {
    match base.to_ascii_uppercase() {
        '=' => 0,  'A' => 1,  'C' => 2,  'M' => 3,
        'G' => 4,  'R' => 5,  'S' => 6,  'V' => 7,
        'T' => 8,  'W' => 9,  'Y' => 10, 'H' => 11,
        'K' => 12, 'D' => 13, 'B' => 14, 'N' => 15,
        _   => 15,
    }
}


/// Converts BAM file to GBAM file. Sorts BAM file in process. This uses the `bam_parallel` reader.
pub fn bam_sort_to_gbam(in_path: &str, out_path: &str, codec: Codecs, mut sort_temp_mode: Option<String>, temp_dir: Option<PathBuf>, full_command: String, index_sort: bool) {
    let fin_for_ref_seqs = File::open(in_path).expect("failed");
    
    let mut reader_for_header_only = Reader::new(fin_for_ref_seqs, 1, None);
    let (sam_header, ref_seqs, _) =
        read_sam_header_and_ref_seqs(&mut reader_for_header_only);


    let fin = File::open(in_path).expect("failed");
    let fout = File::create(out_path).expect("failed");

    let file_size = fin.metadata().unwrap().len();

    let buf_reader = BufReader::new(fin);
    let buf_writer = BufWriter::new(fout);

    let mut writer = Writer::new(
        buf_writer,
        vec![codec; FIELDS_NUM],
        8,
        vec![Fields::RefID],
        ref_seqs,
        sam_header,
        full_command,
        true
    );

    let tmp_dir_path = temp_dir.map_or(std::env::temp_dir(), |path| path);
    if sort_temp_mode.is_none() {
        sort_temp_mode = Some(String::from_str("file").unwrap());
    }
    let tmp_medium_mode = match sort_temp_mode.unwrap().as_str() {
        "file" => TempFilesMode::RegularFiles,
        "lz4_file" => TempFilesMode::LZ4CompressedFiles,
        "ram" => TempFilesMode::InMemoryBlocks,
        "lz4_ram" => TempFilesMode::InMemoryBlocksLZ4,
        _ => panic!("Unknown sort_temp_mode mode."),
    };
    
    let index_file = if index_sort {
        Some(BufWriter::with_capacity(33_554_432, File::create(out_path.clone().to_owned()+".gbai").unwrap()))
    }
    else{None};

    let dir = TempDir::new_in(tmp_dir_path, "BAM sort temporary directory.").unwrap();

    sort::sort_bam(
        MEM_LIMIT,
        buf_reader,
        &mut writer,
        &dir,
        0,
        8,
        tmp_medium_mode,
        index_file,
        sort::SortBy::CoordinatesAndStrand,
        Some(file_size)
    )
    .unwrap();

    writer.finish().unwrap();
}

/// Consumes SAM header from input BAM reader.
///
///
/// # Returns
/// Tuple of 3 elements.
///
/// **tuple.0** -> all bytes of BAM header including reference sequences.
///
/// **tuple.1** -> parsed reference sequences from BAM header.
///
/// **tuple.2** -> offset to reference sequences in tuple.0. It's before n_ref uint32_t.
fn read_sam_header_and_ref_seqs(reader: &mut Reader) -> (Vec<u8>, Vec<(String, u32)>, usize) {
    let (bytes_of_header, ref_sequences_offset) = reader.read_header().unwrap();
    let sequences = parse_reference_sequences(&bytes_of_header[ref_sequences_offset..]).unwrap();
    (bytes_of_header, sequences, ref_sequences_offset)
}

fn get_bam_reader_gbam_writer(
    in_path: &str,
    out_path: &str,
    codec: Codecs,
    full_command: String,
) -> (Reader, Writer<BufWriter<File>>) {
    let fin = File::open(in_path).expect("failed");
    let fout = File::create(out_path).expect("failed");

    let file_size = fin.metadata().unwrap().len();

    let buf_reader = BufReader::new(fin);
    let buf_writer = BufWriter::new(fout);

    let mut bgzf_reader = Reader::new(buf_reader, 4, Some(file_size));

    let (sam_header, ref_seqs, _) = read_sam_header_and_ref_seqs(&mut bgzf_reader);

    let writer = Writer::new(
        buf_writer,
        vec![codec; FIELDS_NUM],
        8,
        vec![Fields::RefID],
        ref_seqs,
        sam_header,
        full_command,
        false,
    );

    (bgzf_reader, writer)
}

fn get_sam_reader_gbam_writer1(
    in_path: &str,
    out_path: &str,
    codec: Codecs,
    full_command: String,
) -> (SamReader<BufReader<File>>, Writer<BufWriter<File>>) {
    let fin = File::open(in_path).expect("failed");
    let fout = File::create(out_path).expect("failed");

    let file_size = fin.metadata().unwrap().len();

    let buf_reader = BufReader::new(fin);
    let buf_writer = BufWriter::new(fout);

    // Create the SAM reader (parses header inside `new`)
    let sam_reader = SamReader::new(buf_reader);

    // Extract SAM header bytes and ref seqs
    let (sam_header, _) = sam_reader.header().expect("failed to get header");
    let ref_seqs = sam_reader.reference_sequences();

    // Create GBAM writer
    let writer = Writer::new(
        buf_writer,
        vec![codec; FIELDS_NUM],
        8,
        vec![Fields::RefID],
        ref_seqs,
        sam_header,
        full_command,
        false,
    );

    (sam_reader, writer)
}

fn get_sam_reader_gbam_writer(
    in_path: Option<&str>, // now optional
    out_path: &str,
    codec: Codecs,
    full_command: String,
) -> (SamReader<BufReader<Box<dyn io::Read>>>, Writer<BufWriter<File>>) {
    // Input: file or stdin
    let reader: Box<dyn Read> = match in_path {
        Some(path) => Box::new(File::open(path).expect("failed to open input file")),
        None => Box::new(io::stdin()), // <-- default to stdin
    };
    let buf_reader = BufReader::new(reader);

    // Output: file
    let fout = File::create(out_path).expect("failed to create output file");
    let buf_writer = BufWriter::new(fout);

    // Create the SAM reader
    let sam_reader = SamReader::new(buf_reader);
    let (sam_header, _) = sam_reader.header().expect("failed to get header");
    let ref_seqs = sam_reader.reference_sequences();

    let writer = Writer::new(
        buf_writer,
        vec![codec; FIELDS_NUM],
        8,
        vec![Fields::RefID],
        ref_seqs,
        sam_header,
        full_command,
        false,
    );

    (sam_reader, writer)
}
