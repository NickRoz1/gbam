use super::{
    compression_enum_size, u8_size, Compression, Fields, RawRecord, FIELDS_NUM, GBAM_MAGIC,
};
use byteorder::{LittleEndian, ReadBytesExt, WriteBytesExt};
use flate2::write::ZlibEncoder;
use std::io;
use std::io::Write;
use std::path::Path;
use zstd::stream::Encoder;

/// A GBAM writer.
pub struct Writer<W: Write> {
    writer: W,
    cur_block_offset: u64,
}

// To leverage column oriented storage the data (iterator) can be requested as
// special tuple and only the necessary columns will be traversed.
//
// The rowgroups are of variable size. If variable sized field overflow maximum
// rowgroup size, it placed in a queue for later processing, and rowgroup
// without this record is processed.

impl<W> Writer<W>
where
    W: Write,
{
    /// Writes ID data to passed writers, wraps them into encoder and creates
    /// Writer object using them.
    pub fn new(mut inner: W) -> io::Result<Self> {
        inner.write(&GBAM_MAGIC)?;
        Ok(Self {
            writer: inner,
            cur_block_offset: 0,
        })
    }
}

impl<W> Writer<W>
where
    W: Write,
{
    /// Splits the record field by field and writes them to separate files.
    /// Works correctly only when full record is provided.
    pub fn write_record(&mut self, rec: &RawRecord) -> io::Result<()> {
        for field in Fields::iterator() {
            match field {
                Fields::LName | Fields::SequenceLength | Fields::NCigar | Fields::RawTagsLen => {
                    let offset = &mut self.offsets[*field as usize];
                    match field {
                        Fields::LName => self.encoded_writers[*field as usize]
                            .write_u32::<LittleEndian>(*offset as u32)?,
                        Fields::SequenceLength => self.encoded_writers[*field as usize]
                            .write_u32::<LittleEndian>(*offset as u32)?,
                        Fields::NCigar => self.encoded_writers[*field as usize]
                            .write_u32::<LittleEndian>(*offset as u32)?,
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
