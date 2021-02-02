use super::{Compression, Fields, RawRecord, FIELDS_NUM, GBAM_MAGIC};
use byteorder::{LittleEndian, WriteBytesExt};
use flate2::write::ZlibEncoder;
use std::io;
use std::io::Write;
use std::path::Path;
use zstd::stream::Encoder;

/// A GBAM writer.
pub struct Writer<'a> {
    // One writer per record field type.
    encoded_writers: Vec<Box<dyn Write + 'a>>,
    offsets: Vec<u64>,
}

impl<'a> Writer<'a> {
    /// Writes ID data to passed writers, wraps them into encoder and creates
    /// Writer object using them.
    pub fn new<W: Write + 'a>(mut inners: Vec<W>, compr_type: Compression) -> io::Result<Self> {
        Self::set_file_meta(&mut inners)?;
        let mut encoded_writers: Vec<Box<dyn Write + 'a>> = Vec::new();
        for writer in inners {
            match compr_type {
                Compression::ZSTD => {
                    encoded_writers.push(Box::new(Encoder::new(writer, 1)?.auto_finish()))
                }
                Compression::FLATE2 => encoded_writers.push(Box::new(ZlibEncoder::new(
                    writer,
                    flate2::Compression::default(),
                ))),
            }
        }
        Ok(Self {
            encoded_writers: encoded_writers,
            offsets: vec![0; FIELDS_NUM as usize],
        })
    }

    /// Writes magic number in each writer and follows with record field number.
    fn set_file_meta<W: Write + 'a>(inners: &mut Vec<W>) -> io::Result<()> {
        if inners.len() != FIELDS_NUM {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "There has to be one writer per each record field.",
            ));
        }
        for (inner, field_type) in inners.iter_mut().zip(Fields::iterator()) {
            inner.write(&GBAM_MAGIC)?;
            inner.write(&[*field_type as u8])?;
            inner.write(b"\x01")?;
        }
        Ok(())
    }
}

impl<'a> Writer<'a> {
    /// Splits the record field by field and writes them to separate files.
    /// Works correctly only when full record is provided.
    pub fn write_record(&mut self, rec: &RawRecord) -> io::Result<()> {
        for field in Fields::iterator() {
            match field {
                Fields::LName | Fields::SequenceLength | Fields::NCigar | Fields::RawTagsLen => {
                    let offset = &mut self.offsets[*field as usize];
                    match field {
                        Fields::LName => {
                            self.encoded_writers[*field as usize].write_u8(*offset as u8)?
                        }
                        Fields::SequenceLength => self.encoded_writers[*field as usize]
                            .write_u32::<LittleEndian>(*offset as u32)?,
                        Fields::NCigar => self.encoded_writers[*field as usize]
                            .write_u16::<LittleEndian>(*offset as u16)?,
                        Fields::RawTagsLen => self.encoded_writers[*field as usize]
                            .write_u32::<LittleEndian>(*offset as u32)?,
                        _ => panic!("This field is not supported: {} \n", *field as usize),
                    }
                    *offset += rec.get_len_val(field) as u64;
                }
                _ => {
                    self.encoded_writers[*field as usize].write(rec.get_bytes(field))?;
                }
            }
        }
        Ok(())
    }
    // /// Flush the buffers.
    // pub fn flush(&mut self) -> io::Result<()> {}
}

impl<'a> Drop for Writer<'a> {
    fn drop(&mut self) {
        for writer in &mut self.encoded_writers {
            writer.flush().expect("Closing Writer failed.");
        }
    }
}
