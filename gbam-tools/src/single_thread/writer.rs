use super::meta::{BlockMeta, FileInfo, FileMeta};
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

impl<W> Writer<W>
where
    W: Write + Seek,
{
    pub fn new(mut inner: W) -> Self {
        // Make space for the FileInfo to be written into.
        inner
            .seek(SeekFrom::Start((u64_size + u32_size * 2 + u64_size) as u64))
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
            if *offset > SIZE_LIMIT {
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
        let mut seek_pos = 0;
        if !field_meta.is_empty() {
            let previous = field_meta.last().unwrap();
            seek_pos = previous.seekpos + (previous.numitems * item_size) as u64;
        }
        // Space taken by values
        let offset = self.offsets[*field as usize];
        BlockMeta {
            seekpos: seek_pos,
            numitems: self.num_items[*field as usize],
        }
    }

    /// Terminates the writer. Always call after writting all the data.
    pub fn finish(&mut self) {
        // DONT DELETE THIS! THIS HOW IT SUPPOSED TO WORK WHEN ALL FIELDS ARE AVAILABLE!
        // for field in self.fields_to_flush.iter_mut() {
        //     *field = true;
        // }
        self.fields_to_flush[Fields::Mapq as usize] = true;
        self.fields_to_flush[Fields::Pos as usize] = true;
        self.flush();
        let meta_start_pos = self.inner.seek(SeekFrom::Current(0)).unwrap();
        // Write meta
        let main_meta = serde_json::to_string(&self.file_meta).unwrap();
        self.inner.write(&main_meta.as_bytes()[..]);
        // Revert back to the beginning of the file
        self.inner.seek(SeekFrom::Start(0)).unwrap();
        let file_meta = FileInfo::new([1, 0], meta_start_pos);
        self.inner.write(&Into::<Vec<u8>>::into(file_meta)[..]);
    }
}
