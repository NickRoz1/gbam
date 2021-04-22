use super::meta::{BlockMeta, FileInfo, FileMeta, FILE_INFO_SIZE};
use super::SIZE_LIMIT;
use crate::{u32_size, u64_size, u8_size, Fields, RawRecord, FIELDS_NUM};
use byteorder::{LittleEndian, WriteBytesExt};
use std::io::{Seek, SeekFrom, Write};
static GBAM_MAGIC: &[u8] = b"geeBAM10";

/// Groups records before writing out to file.
pub struct Writer<W>
where
    W: Write + Seek,
{
    chunks: Vec<Vec<u8>>,
    // Current item index
    offsets: [usize; FIELDS_NUM],
    num_items: [u32; FIELDS_NUM],
    fields_to_flush: [bool; FIELDS_NUM],
    file_meta: FileMeta,
    inner: W,
}
static mut num: u32 = 0;
impl<W> Writer<W>
where
    W: Write + Seek,
{
    pub fn new(mut inner: W) -> Self {
        // Make space for the FileInfo to be written into.
        inner
            .seek(SeekFrom::Start((FILE_INFO_SIZE) as u64))
            .unwrap();
        Writer {
            chunks: vec![vec![0; SIZE_LIMIT]; FIELDS_NUM],
            offsets: [0; FIELDS_NUM],
            num_items: [0; FIELDS_NUM],
            fields_to_flush: [false; FIELDS_NUM],
            file_meta: FileMeta::new(),
            inner: inner,
        }
    }
    pub fn push_record(&mut self, record: &RawRecord) {
        for field in Fields::iterator().filter(|v| **v == Fields::Mapq || **v == Fields::Pos) {
            let cur_chunk = &mut self.chunks[*field as usize];
            let offset = &mut self.offsets[*field as usize];
            let item_counter = &mut self.num_items[*field as usize];
            let new_data = record.get_bytes(field);
            cur_chunk[*offset..*offset + new_data.len()].clone_from_slice(new_data);
            *offset += new_data.len();
            *item_counter += 1;
        }

        let mut flush_required = false;
        for (idx, offset) in self.offsets.iter_mut().enumerate() {
            if *offset >= SIZE_LIMIT {
                self.fields_to_flush[idx] = true;
                flush_required = true;
            }
        }

        if flush_required {
            self.flush()
        }
    }

    fn flush(&mut self) {
        let fields_to_flush = self.fields_to_flush.clone();
        for (_, field) in fields_to_flush
            .iter()
            .zip(Fields::iterator())
            .filter(|(v, field)| **v)
        {
            let meta = self.generate_meta(field);
            let field_meta = self.file_meta.get_blocks(field);
            field_meta.push(meta);
            // Write the data
            self.inner
                .write(&self.chunks[*field as usize][0..self.offsets[*field as usize]]);
        }
        for (idx, field) in self.fields_to_flush.iter_mut().enumerate() {
            if *field {
                self.offsets[idx] = 0;
                self.num_items[idx] = 0;
            }
            *field = false;
        }
    }

    fn generate_meta(&mut self, field: &Fields) -> BlockMeta {
        let item_size = self.file_meta.get_field_size(field);
        let field_meta = self.file_meta.get_blocks(field);
        let mut seek_pos = self.inner.seek(SeekFrom::Current(0 as i64)).unwrap();
        // Space taken by values
        let offset = self.offsets[*field as usize];
        BlockMeta {
            seekpos: seek_pos,
            numitems: self.num_items[*field as usize],
        }
    }

    /// Terminates the writer. Always call after writting all the data. Returns
    /// total amount of bytes written.
    pub fn finish(&mut self) -> std::io::Result<u64> {
        // DONT DELETE THIS! THIS HOW IT SUPPOSED TO WORK WHEN ALL FIELDS ARE AVAILABLE!
        // for field in self.fields_to_flush.iter_mut() {
        //     *field = true;
        // }
        self.fields_to_flush[Fields::Mapq as usize] = true;
        self.fields_to_flush[Fields::Pos as usize] = true;
        self.flush();
        let meta_start_pos = self.inner.seek(SeekFrom::Current(0))?;
        // Write meta
        let main_meta = serde_json::to_string(&self.file_meta).unwrap();
        self.inner.write(&main_meta.as_bytes()[..])?;

        let total_bytes_written = self.inner.seek(SeekFrom::Current(0))?;
        // Revert back to the beginning of the file
        self.inner.seek(SeekFrom::Start(0)).unwrap();
        let file_meta = FileInfo::new([1, 0], meta_start_pos);
        self.inner.write(&Into::<Vec<u8>>::into(file_meta)[..])?;
        Ok(total_bytes_written)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::single_thread::reader::*;
    use byteorder::ReadBytesExt;
    use std::io::Cursor;
    #[test]
    fn test_writer() {
        let raw_records = vec![RawRecord::default(); 2];
        let mut buf: Vec<u8> = vec![0; SIZE_LIMIT];
        let out = Cursor::new(&mut buf[..]);
        let mut writer = Writer::new(out);
        for rec in raw_records.iter() {
            writer.push_record(rec);
        }
        let total_bytes_written = writer.finish().unwrap();
        buf.resize(total_bytes_written as usize, 0);

        let in_cursor = Box::new(Cursor::new(buf));
        let mut parsing_template = ParsingTemplate::new();
        parsing_template.set_all();
        let mut reader = Reader::new(in_cursor, parsing_template).unwrap();
        let mut it = raw_records.iter();
        while let Some(rec) = reader.next() {
            let rec_orig = it.next().unwrap();
            let orig_map_q = rec_orig.get_bytes(&Fields::Mapq)[0];
            let orig_pos = rec_orig
                .get_bytes(&Fields::Pos)
                .read_u32::<LittleEndian>()
                .unwrap();
            assert_eq!(rec.pos.unwrap(), orig_pos);
            assert_eq!(rec.mapq.unwrap(), orig_map_q);
        }
    }
}
