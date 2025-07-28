use std::convert::{TryFrom, TryInto};
use std::fs::File;
use std::io::BufWriter;
use std::num::NonZeroUsize;

use noodles_core::position::Position;
use noodles_bam::lazy::record::Data as LazyData;
use noodles_bam::writer::Writer as BamWriter;
use noodles_sam::header::{self, record::value::{map::ReferenceSequence, Map}};
use noodles_sam::record::{
    ReferenceSequenceName, Sequence, QualityScores, ReadName, MappingQuality, Data,
    data::field::{Tag, Value, Type, value::{Array, Character}}
};

use noodles_sam::alignment::{Record as SamRecord, record::Builder};

use crate::reader::parse_tmplt::ParsingTemplate;
use crate::reader::reader::Reader;
use crate::reader::records::Records;

/// Converts GBAM file to BAM file. This uses the `noodles bam writer`.
pub fn gbam_to_bam(in_path: &str, out_path: &str) -> Result<(), Box<dyn std::error::Error>> {
    let file = File::open(in_path).unwrap();
    let mut template = ParsingTemplate::new();
    template.set_all();
    let mut reader = Reader::new(file, template).unwrap();

    // Build SAM header
    let mut header_builder = header::Header::builder()
        .set_header(Default::default());
    let ref_seqs = reader.file_meta.get_ref_seqs();

    // Contruct the header
    for (name, len) in ref_seqs {
        let name: ReferenceSequenceName = name.parse()?; // convert String to refseq name
        let len = NonZeroUsize::new(*len as usize).unwrap(); // safe if len > 0
        let mut map = Map::<ReferenceSequence>::new(len);

        header_builder = header_builder.add_reference_sequence(name, map);
    }

    let header = header_builder.build();

    // Create output BAM file and write the header
    let out_file = File::create(out_path)?;
    let mut writer = BamWriter::new(BufWriter::new(out_file));

    writer.write_header(&header)?;

    let mut records_it = Records::new(&mut reader);
    while let Some(gbam_record) = records_it.next_rec() {
        // Use the Builder to create a SAM record
        let mut builder = Builder::default();

        // Set the read name
        if let Some(read_name) = &gbam_record.read_name {
            // Remove null byte (if present) from the read name
            let sanitized_read_name = read_name
                .iter()
                .copied()
                .take_while(|&b| b != 0)
                .collect::<Vec<u8>>();
        
            if let Ok(read_name) = ReadName::try_from(sanitized_read_name) {
                builder = builder.set_read_name(read_name);
            } else {
                eprintln!("Invalid read name: {:?}", gbam_record.read_name);
            }
        }

        // Set the reference sequence ID
        if let Some(refid) = gbam_record.refid {
            if refid != -1 {
                builder = builder.set_reference_sequence_id(refid as usize);
            }
        }

        // Set the position
        if let Some(mut pos) = gbam_record.pos {
            // Add 1 to pos to match the original BAM
            pos += 1;
        
            // Ensure the position is valid before setting it
            if let Some(position) = Position::new(pos as usize) {
                builder = builder.set_alignment_start(position);
            } else {
                eprintln!("Invalid position: {}", pos);
                // Optionally handle the invalid case (e.g., skip setting the alignment start)
            }
        }

        // Set the Fags
        builder = builder
            .set_flags(gbam_record.flag.unwrap().into());

        // Set the mapping quality
        if let Some(mapq) = MappingQuality::new(gbam_record.mapq.unwrap()) {
            builder = builder.set_mapping_quality(mapq);
        }

        // Set the template length
        builder = builder.set_template_length(gbam_record.tlen.unwrap() as i32);

        // Set the next ref ID
        if let Some(next_ref_id) = gbam_record.next_ref_id {
            if next_ref_id != -1 {
                builder = builder.set_mate_reference_sequence_id(next_ref_id as usize);
            }
        }

        // Set the next position
        if let Some(next_pos) = gbam_record.next_pos {
            if next_pos != -1 {
                let mut adjusted_next_pos = next_pos + 1; // Create a mutable copy and add 1
                if let Some(position) = Position::new(adjusted_next_pos as usize) {
                    builder = builder.set_mate_alignment_start(position);
                } else {
                    eprintln!("Invalid next_pos: {}", adjusted_next_pos);
                    // Handle the invalid case (e.g., skip setting the mate alignment start)
                }
            }
        }

        // Set the sequence
        if let Some(seq) = &gbam_record.seq {
            let seq_bytes = seq.clone().into_bytes();
            if let Ok(sequence) = Sequence::try_from(seq_bytes) {
                builder = builder.set_sequence(sequence);
            } 
        }

        // Set the quality scores
        if let Some(qual) = &gbam_record.qual {
            // Convert Vec<u8> to Vec<Score>
            let scores: Vec<noodles_sam::record::quality_scores::Score> = qual
                .iter()
                .map(|&q| noodles_sam::record::quality_scores::Score::try_from(q).unwrap())
                .collect();

            // Create QualityScores from Vec<Score>
            let quality_scores = QualityScores::from(scores);

            builder = builder.set_quality_scores(quality_scores);
        }

        // Set the cigar string
        if let Some(cigar) = &gbam_record.cigar {
            let converted_cigar: noodles_sam::record::Cigar = cigar
                .ops()
                .map(|op| {
                    // Convert query::cigar::Op to noodles_sam::record::cigar::Op
                    let kind = match op.op_type() {
                        'M' => noodles_sam::record::cigar::op::Kind::Match,
                        'I' => noodles_sam::record::cigar::op::Kind::Insertion,
                        'D' => noodles_sam::record::cigar::op::Kind::Deletion,
                        'N' => noodles_sam::record::cigar::op::Kind::Skip,
                        'S' => noodles_sam::record::cigar::op::Kind::SoftClip,
                        'H' => noodles_sam::record::cigar::op::Kind::HardClip,
                        'P' => noodles_sam::record::cigar::op::Kind::Pad,
                        '=' => noodles_sam::record::cigar::op::Kind::SequenceMatch,
                        'X' => noodles_sam::record::cigar::op::Kind::SequenceMismatch,
                        _ => panic!("Unexpected CIGAR operation type: {}", op.op_type()),
                    };
                    noodles_sam::record::cigar::Op::new(kind, op.length().try_into().unwrap())
                })
                .collect();

            builder = builder.set_cigar(converted_cigar);
        }

        // Set the optional tags
        if let Some(tags) = &gbam_record.tags {
            // Parse the raw binary tags into noodles_sam::record::Data
            match parse_tags(tags) {
                Ok(data) => {
                    builder = builder.set_data(data);
                }
                Err(e) => {
                    eprintln!("Error parsing tags: {}", e);
                }
            }
        }

        let mut sam_record = builder.build();

        if let Err(e) = writer.write_record(&header, &sam_record) {
            // Capture the backtrace
            let backtrace = std::backtrace::Backtrace::capture();
        
            // Check the backtrace status
            if backtrace.status() == std::backtrace::BacktraceStatus::Captured {
                eprintln!("Backtrace:\n{:?}", backtrace);
            }
        
            return Err(e.into());
        }
    }

    Ok(())
}

fn parse_tags(tags: &[u8]) -> Result<Data, Box<dyn std::error::Error>> {
    let mut data = Data::default();
    let mut i = 0;

    while i < tags.len() {
        // Parse the tag (2 bytes)
        let tag = Tag::try_from(&tags[i..i + 2])?;
        i += 2;

        // Parse the type (1 byte)
        let ty = tags[i] as char;
        i += 1;

        // Parse the value based on the type
        let value = match ty {
            'A' => {
                // Character
                let c = tags[i] as char;
                i += 1;
                Value::Character(Character::try_from(c)?)
            }
            // For a tag which has [78, 77, 67, 1] bytes, it is decoded as NM:C:1.
            // But here type C is not a valid BAM tag type. Hence if the tag type is C
            // value should be converted to 8-bit signed integer type.
            'C' => {
                // Integer
                let val = u8::from_le_bytes(tags[i..i + 1].try_into()?);
                i += 1;
                Value::UInt8(val)
            }
            'c' => {
                // Integer
                let val = i8::from_le_bytes(tags[i..i + 1].try_into()?);
                i += 1;
                Value::Int8(val)
            }
            'i' => {
                // Integer
                let val = i32::from_le_bytes(tags[i..i + 4].try_into()?);
                i += 4;
                Value::Int32(val)
            }
            'f' => {
                // Float
                let val = f32::from_le_bytes(tags[i..i + 4].try_into()?);
                i += 4;
                Value::Float(val)
            }
            'Z' => {
                // String (null-terminated)
                let end = tags[i..].iter().position(|&b| b == 0).ok_or("Missing null terminator")?;
                let s = std::str::from_utf8(&tags[i..i + end])?.to_string();
                i += end + 1;
                Value::String(s)
            }
            'B' => {
                // Array
                let sub_type = tags[i] as char;
                i += 1;

                let count = u32::from_le_bytes(tags[i..i + 4].try_into()?);
                i += 4;

                let mut values = Vec::new();

                // Push the ASCII value of the sub_type into the values array
                values.push(sub_type as u8);

                for _ in 0..count {
                    match sub_type {
                        'C' => {
                            let val = tags[i];
                            i += 1;
                            values.push(val);
                        }
                        'c' => {
                            let val = tags[i] as i8;
                            i += 1;
                            values.push(val as u8); // Convert to u8 for consistency
                        }
                        'S' => {
                            let val = u16::from_le_bytes(tags[i..i + 2].try_into()?);
                            i += 2;
                            values.extend_from_slice(&val.to_le_bytes());
                        }
                        's' => {
                            let val = i16::from_le_bytes(tags[i..i + 2].try_into()?);
                            i += 2;
                            values.extend_from_slice(&val.to_le_bytes());
                        }
                        'I' => {
                            let val = u32::from_le_bytes(tags[i..i + 4].try_into()?);
                            i += 4;
                            values.extend_from_slice(&val.to_le_bytes());
                        }
                        'i' => {
                            let val = i32::from_le_bytes(tags[i..i + 4].try_into()?);
                            i += 4;
                            values.extend_from_slice(&val.to_le_bytes());
                        }
                        'f' => {
                            let val = f32::from_le_bytes(tags[i..i + 4].try_into()?);
                            i += 4;
                            values.extend_from_slice(&val.to_le_bytes());
                        }
                        _ => return Err(format!("Unsupported B array subtype: {}", sub_type).into()),
                    }
                }

                // Construct the appropriate Array variant
                let array = match sub_type {
                    'C' => Array::UInt8(values),
                    'c' => Array::Int8(values.into_iter().map(|v| v as i8).collect()),
                    'S' => Array::UInt16(values.chunks(2).map(|chunk| u16::from_le_bytes(chunk.try_into().unwrap())).collect()),
                    's' => Array::Int16(values.chunks(2).map(|chunk| i16::from_le_bytes(chunk.try_into().unwrap())).collect()),
                    'I' => Array::UInt32(values.chunks(4).map(|chunk| u32::from_le_bytes(chunk.try_into().unwrap())).collect()),
                    'i' => Array::Int32(values.chunks(4).map(|chunk| i32::from_le_bytes(chunk.try_into().unwrap())).collect()),
                    'f' => Array::Float(values.chunks(4).map(|chunk| f32::from_le_bytes(chunk.try_into().unwrap())).collect()),
                    _ => return Err(format!("Unsupported B array subtype: {}", sub_type).into()),
                };

                Value::Array(array)
            }
            _ => return Err(format!("Unsupported tag type: {}", ty).into()),
        };

        // Insert the tag-value pair into the Data object
        data.insert(tag, value);
    }

    Ok(data)
}

// pub fn gbam_to_bam(in_path: &str, out_path: &str) {
//     let file = File::open(in_path).unwrap();
//     let mut template = ParsingTemplate::new();
//     template.set_all();
//     let mut reader = crate::reader::reader::Reader::new(file, template).unwrap();

//     let mut bam_header = bam::Header::new();
//     let ref_seqs = reader.file_meta.get_ref_seqs();

//     for ref_seq in ref_seqs {
//         bam_header.push_record(
//             bam::header::HeaderRecord::new(b"SQ")
//                 .push_tag(b"SN", &ref_seq.0)
//                 .push_tag(b"LN", &ref_seq.1),
//         );
//     }

//     let mut records_it = Records::new(&mut reader);

//     let mut out = bam::Writer::from_path(out_path, &bam_header, bam::Format::Bam).unwrap();
//     out.set_threads(4).unwrap();

//     let mut cigar_buf = Vec::new();
//     let mut record_count = 0;
//     while let Some(rec) = records_it.next_rec() {
//         let mut record = bam::Record::new();
//         record_count += 1;
//         if record_count == 1 {
//             println!("Raw tags: {:?}", rec.tags);
//         }

//         record.set_bin(rec.bin.unwrap());
//         record.set_tid(rec.refid.unwrap());
//         record.set_mapq(rec.mapq.unwrap());
//         record.set_pos(rec.pos.unwrap() as i64);
//         record.set_flags(rec.flag.unwrap());
//         record.set_mtid(rec.next_ref_id.unwrap());
//         record.set_mpos(rec.next_pos.unwrap() as i64);
//         record.set_insert_size(rec.tlen.unwrap() as i64);
//         let rec_seq_len = rec.seq.as_ref().unwrap().len();
//         let mut qual = rec.qual.as_ref().unwrap().clone();
//         if qual.is_empty() {
//             qual = vec![255; rec_seq_len];
//         }

//         cigar_buf.clear();
//         rec.cigar.as_ref().unwrap().ops().for_each(|op| {
//             cigar_buf
//                 .write_all(op.length().to_string().as_bytes())
//                 .unwrap();
//             cigar_buf.push(op.op_type() as u8);
//         });

//         let bam_cigar = bam::record::CigarString::try_from(&cigar_buf[..]).unwrap();
//         record.set_data(&rec.tags.as_ref().unwrap()[..]);
//         record.set(
//             &rec.read_name.as_ref().unwrap()[..rec.read_name.as_ref().unwrap().len() - 1],
//             Some(&bam_cigar),
//             rec.seq.as_ref().unwrap().as_bytes(),
//             &qual[..],
//         );

//         out.write(&record).unwrap();
//     }
// }
