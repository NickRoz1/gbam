use super::meta::{BlockMeta, Codecs, FileInfo, FileMeta, FILE_INFO_SIZE, Stat};
use crate::compressor::{CompressTask, Compressor, OrderingKey};
use crate::{SIZE_LIMIT, U32_SIZE};
use bam_tools::record::bamrawrecord::BAMRawRecord;
use bam_tools::record::fields::{
    field_type, is_data_field, var_size_field_to_index, FieldType, Fields, FIELDS_NUM,
};
use byteorder::{LittleEndian, WriteBytesExt, ReadBytesExt};
use crc32fast::Hasher;
use std::borrow::Cow;
use std::convert::TryInto;
use std::convert::TryFrom;
use std::io::{Seek, SeekFrom, Write};
use std::collections::HashMap;
use once_cell::sync::Lazy;
use std::fs;
use std::path::Path;
use serde_json::Value;
use std::str::FromStr;

pub(crate) struct BlockInfo {
    pub numitems: u32,
    pub uncompr_size: usize,
    pub field: Fields,
    // Interpretation is up to the reader.
    pub stats: Option<Stat>,
    pub codec: Codecs,
}

impl Default for BlockInfo {
    fn default() -> Self {
        Self {
            numitems: 0,
            uncompr_size: 0,
            field: Fields::RefID,
            stats: None,
            codec: Codecs::Brotli,
        }
    }
}

pub static FIELD_CODEC_MAP: Lazy<HashMap<Fields, Codecs>> = Lazy::new(|| {
    let path = Path::new("codec_map.json");
    let json_str = fs::read_to_string(path).expect("Failed to read codec_map.json");

    let parsed: Value = serde_json::from_str(&json_str).expect("Invalid JSON format");

    let mut map = HashMap::new();
    for (field_str, codec_str) in parsed.as_object().unwrap() {
        let field = Fields::from_str(field_str)
            .unwrap_or_else(|_| panic!("Unknown field in JSON: {}", field_str));
        let codec = match codec_str.as_str().unwrap() {
            "Brotli" => Codecs::Brotli,
            "Zstd" => Codecs::Zstd,
            "Lz4" => Codecs::Lz4,
            "Gzip" => Codecs::Gzip,
            "NoCompression" => Codecs::NoCompression,
            other => panic!("Unsupported codec: {}", other),
        };
        map.insert(field, codec);
    }
    map
});

/// The data is held in blocks.
///
/// Fixed sized fields are written as fixed size blocks into file. All blocks
/// (for fixed size fields) except last one contain equal amount of data.
///
/// Variable sized fields are written as fixed size blocks. Blocks may contain
/// different amount of data. Variable sized fields are accompanied by separate
/// index in separate block for fixed size fields. Groups records before writing
/// out to file.
pub struct Writer<WS>
where
    WS: Write + Seek,
{
    file_info: FileInfo,
    file_meta: FileMeta,
    columns: Vec<Box<dyn Column>>,
    compressor: Compressor,
    inner: WS,
}

impl<WS> Writer<WS>
where
    WS: Write + Seek,
{
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        mut inner: WS,
        codecs: Vec<Codecs>,
        thread_num: usize,
        collect_stats_for: Vec<Fields>,
        ref_seqs: Vec<(String, u32)>,
        sam_header: Vec<u8>,
        full_command: String,
        is_sorted: bool,
    ) -> Self {
        inner
            .seek(SeekFrom::Start((FILE_INFO_SIZE) as u64))
            .unwrap();

        let mut columns = Vec::new();

        let mut count = 0;
        for field in Fields::iterator().filter(|f| is_data_field(f)) {
            let stat_collector = collect_stats_for.iter().find(|f| *f == field).and(Some(Stat::default()));
            let col = match field_type(field) {
                FieldType::FixedSized => {
                    Box::new(FixedColumn::new(*field, stat_collector)) as Box<dyn Column>
                }
                FieldType::VariableSized => {
                    // Index column +1.
                    count += 1;
                    Box::new(VariableColumn::new(*field, stat_collector)) as Box<dyn Column>
                }
            };
            columns.push(col);
            count += 1;
        }
        debug_assert!(count == FIELDS_NUM);

        Self {
            // TODO: Codecs (currently only one is supported).
            file_meta: FileMeta::new(codecs[0], ref_seqs, sam_header),
            inner,
            compressor: Compressor::new(thread_num),
            columns,
            file_info: FileInfo::new([1, 0], 0, 0, full_command, is_sorted),
        }
    }

    pub fn new_no_stats(
        inner: WS,
        codecs: Vec<Codecs>,
        thread_num: usize,
        ref_seqs: Vec<(String, u32)>,
        sam_header: Vec<u8>,
        full_command: String,
        is_sorted: bool,
    ) -> Self {
        Self::new(
            inner,
            codecs,
            thread_num,
            Vec::new(),
            ref_seqs,
            sam_header,
            full_command,
            is_sorted
        )
    }

    /// Push BAM record into this writer
    pub fn push_record(&mut self, record: &BAMRawRecord) {
        // Index fields are not written on their own. They hold index data for variable sized fields.
        for col in self.columns.iter_mut() {
            // Attempt to write data in this column. If the column is full it
            // will return bytes for flushing. While loop is here because
            // variable sized columns also have index columns (fixed size)
            // inside and they might also come full and request flushing
            // simultaneously with containing variable sized field column.
            while let WriteStatus::Full(inner) = col.write_record_field(record) {
                flush_field_buffer(
                    &mut self.inner,
                    &mut self.file_meta,
                    &mut self.compressor,
                    inner,
                );
            }
        }
    }

    /// Terminates the writer. Always call after writting all the data. Returns
    /// total amount of bytes written.
    pub fn finish(&mut self) -> std::io::Result<u64> {
        // Flush leftovers
        let mut columns: Vec<Box<dyn Column>> = self.columns.drain(..).collect();
        for (inner, idx) in columns.iter_mut().map(|col| col.get_inners()) {
            let writer = &mut self.inner;
            let meta = &mut self.file_meta;
            let compress = &mut self.compressor;

            flush_field_buffer(writer, meta, compress, inner);
            if let Some(idx_inner) = idx {
                flush_field_buffer(writer, meta, compress, idx_inner);
            }
        }

        for mut task in self.compressor.finish() {
            if let OrderingKey::Key(key) = task.ordering_key {
                write_data_and_update_meta(&mut self.inner, &mut self.file_meta, key, &mut task);
            }
        }

        let meta_start_pos = self.inner.stream_position()?;
        // Write meta
        let main_meta = serde_json::to_string(&self.file_meta).unwrap();
        let main_meta_bytes = main_meta.as_bytes();
        let crc32 = calc_crc_for_meta_bytes(main_meta_bytes);
        self.inner.write_all(main_meta_bytes)?;

        let total_bytes_written = self.inner.stream_position()?;
        // Revert back to the beginning of the file
        self.inner.seek(SeekFrom::Start(0)).unwrap();
        self.inner.write_all(&[0;FILE_INFO_SIZE]).unwrap();
        self.inner.seek(SeekFrom::Start(0)).unwrap();
        let file_info = & mut self.file_info;
        file_info.seekpos = meta_start_pos;
        file_info.crc32 = crc32;
        let file_info_bytes = serde_json::to_string(&file_info).unwrap();
        self.inner.write_all(file_info_bytes.as_bytes())?;
        Ok(total_bytes_written)
    }
}

fn flush_field_buffer<WS: Write + Seek>(
    writer: &mut WS,
    file_meta: &mut FileMeta,
    compressor: &mut Compressor,
    inner: &mut Inner,
) {
    // Use an empty buffer to start the flushing process
    // Don't worry, Vec::new() is temporary, it won't need to fully allocate the Vec as it replaces the reference with the &mut from the reused Buffer
    let data = std::mem::replace(&mut inner.buffer, Vec::new());

    let field = &inner.field;
    let codec = *file_meta.get_field_codec(field);

    compressor.compress_block(
        OrderingKey::Key(inner.block_num),
        inner.generate_block_info(),
        data,
        codec,
    );

    let mut completed_task = compressor.get_compr_block();

    if let OrderingKey::Key(key) = completed_task.ordering_key {
        write_data_and_update_meta(writer, file_meta, key, &mut completed_task);
    }

    // We need to reuse the same buffer for the next task, as it is always the same size so we can avoid re-allocating the same buffer for each processed block
    inner.buffer = completed_task.buf;

    inner.reset_for_new_block();
}

fn write_data_and_update_meta<WS: Write + Seek>(
    writer: &mut WS,
    file_meta: &mut FileMeta,
    key: u64,
    task: &mut CompressTask,
) {
    let compressed_size = task.buf.len();
    let meta = generate_meta(
        writer,
        &mut task.block_info,
        compressed_size.try_into().unwrap(),
    );

    writer.write_all(&task.buf).unwrap();

    let field_meta = file_meta.get_blocks(&task.block_info.field);
    if field_meta.len() <= key as usize {
        field_meta.resize(key as usize + 1, BlockMeta::default());
    }

    // Order as came in
    field_meta[key as usize] = meta;
}

fn generate_meta<S: Seek>(
    writer: &mut S,
    block_info: &mut BlockInfo,
    block_size: u32,
) -> BlockMeta {
    let seekpos = writer.stream_position().unwrap();
    BlockMeta {
        seekpos,
        numitems: block_info.numitems,
        block_size,
        uncompressed_size: block_info.uncompr_size as u64,
        stats: block_info.stats.take(),
    }
}

enum WriteStatus<'a> {
    Written,
    // Column or its index is at capacity. Flush it.
    Full(&'a mut Inner),
}

struct Inner {
    stats_collector: Option<Stat>,
    buffer: Vec<u8>,
    offset: usize,
    field: Fields,
    rec_count: u32,
    block_num: u64,
}

impl Inner {
    pub fn new(field: Fields, stats_collector: Option<Stat>) -> Self {
        Self {
            stats_collector,
            buffer: Vec::new(),
            offset: 0,
            field,
            rec_count: 0,
            block_num: 0,
        }
    }
    pub fn write_data(&mut self, data: &[u8]) -> WriteStatus {
        // At this point everything should be flushed.
        debug_assert!(!self.flush_required(data));

        let limit = std::cmp::max(data.len(), SIZE_LIMIT);
        if self.buffer.len() < limit {
            self.buffer.resize(limit, 0);
        }

        self.buffer[self.offset..self.offset + data.len()].clone_from_slice(data);
        self.offset += data.len();

        self.rec_count += 1;

        WriteStatus::Written
    }

    pub fn flush_required(&self, data: &[u8]) -> bool {
        // At least one record will be written in even if it exceeds SIZE_LIMIT.
        self.offset > 0 && self.offset + data.len() > SIZE_LIMIT
    }

    pub fn reset_for_new_block(&mut self) {
        self.offset = 0;
        self.rec_count = 0;
        self.block_num += 1;
    }

    pub fn generate_block_info(&mut self) -> BlockInfo {
        let stat = if self.stats_collector.is_some(){
            self.stats_collector.replace(Stat::default())
        }
        else{
            None
        };
        let codec = *FIELD_CODEC_MAP.get(&self.field).expect("Missing codec mapping");
        BlockInfo {
            numitems: self.rec_count,
            uncompr_size: self.offset,
            field: self.field,
            stats: stat,
            codec: codec,
        }
    }
}

trait Column {
    // Extracts and writes data from corresponding BAMRawRecord record.
    fn write_record_field(&mut self, rec: &BAMRawRecord) -> WriteStatus;

    fn get_inners(&mut self) -> (&mut Inner, Option<&mut Inner>);
}

/// Column containing fixed sized fields.
struct FixedColumn(Inner);

impl FixedColumn {
    pub fn new(field: Fields, comparator: Option<Stat>) -> Self {
        if comparator.is_some() && field != Fields::RefID && field != Fields::Pos {
            panic!("Stats collection is only supported for RefID and POS fields.");
        }
        Self(Inner::new(field, comparator))
    }
}

impl Column for FixedColumn {
    fn write_record_field(&mut self, rec: &BAMRawRecord) -> WriteStatus {
        let inner = &mut self.0;
        let data = rec.get_bytes(&inner.field);

        if inner.flush_required(data) {
            return WriteStatus::Full(inner);
        }

        if let Some(ref mut stats) = inner.stats_collector {
            stats.update((&data[..]).read_i32::<LittleEndian>().unwrap());
        }

        inner.write_data(data)
    }

    fn get_inners(&mut self) -> (&mut Inner, Option<&mut Inner>) {
        (&mut self.0, None)
    }
}

struct VariableColumn {
    inner: Inner,
    index: FixedColumn,
}

impl VariableColumn {
    pub fn new(field: Fields, comparator: Option<Stat>) -> Self {
        if comparator.is_some() {
            panic!("Stats collection is not supported for variable length fields.");
        }
        Self {
            inner: Inner::new(field, comparator),
            index: FixedColumn::new(var_size_field_to_index(&field), None),
        }
    }
}

impl Column for VariableColumn {
    fn write_record_field(&mut self, rec: &BAMRawRecord) -> WriteStatus {
        let inner = &mut self.inner;
        let index_inner = &mut self.index.0;

        let data = rec.get_bytes(&inner.field);
        let mut idx_buf: [u8; U32_SIZE] = [0; U32_SIZE];

        if index_inner.flush_required(&idx_buf) {
            return WriteStatus::Full(index_inner);
        }

        if inner.flush_required(data) {
            return WriteStatus::Full(inner);
        }

        assert!(inner.stats_collector.is_none());

        inner.write_data(data);
        (&mut idx_buf[..])
            .write_u32::<LittleEndian>(u32::try_from(inner.offset).unwrap())
            .unwrap();
        index_inner.write_data(&idx_buf)
    }

    fn get_inners(&mut self) -> (&mut Inner, Option<&mut Inner>) {
        (&mut self.inner, Some(&mut self.index.0))
    }
}

impl<W> Write for Writer<W>
where
    W: Write + Seek,
{
    /// WARNING: ENSURE THAT BUF CONTAINS A ONE FULL RECORD.
    /// TODO: Implement a proper trait and use it instead of Write in bam_sorting.
    /// Write trait implementation is made to allow passing Write trait objects to sort function in BAM parallel.
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        assert!(!buf.is_empty());
        let wrapper = BAMRawRecord(Cow::Borrowed(buf));
        self.push_record(&wrapper);
        Ok(buf.len())
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

// TODO: Currently end user should manually call finish. Probably can be done
// with a drop. If drop and manual finish used simultaneously, crc32 and meta of
// file will be damaged.
// impl<W> Drop for Writer<W>
// where
//     W: Write + Seek,
// {
//     fn drop(&mut self) {
//         self.finish().unwrap();
//     }
// }

pub(crate) fn calc_crc_for_meta_bytes(bytes: &[u8]) -> u32 {
    let mut hasher = Hasher::new();
    hasher.update(bytes);
    hasher.finalize()
}

// #[ignore]
// #[cfg(test)]
// mod tests {
//     use super::*;
//     use crate::reader::parse_tmplt::*;
//     use crate::reader::reader::*;
//     use byteorder::ReadBytesExt;
//     use std::io::Cursor;
//     #[test]
//     fn test_writer() {
//         // let raw_records = vec![BAMRawRecord::default(); 2];
//         // let mut buf: Vec<u8> = vec![0; SIZE_LIMIT];
//         // let out = Cursor::new(&mut buf[..]);
//         // let mut writer = Writer::new(out, Codecs::Gzip, 8);
//         // for rec in raw_records.iter() {
//         //     writer.push_record(rec);
//         // }
//         // let total_bytes_written = writer.finish().unwrap();
//         // buf.resize(total_bytes_written as usize, 0);

//         // let in_cursor = Box::new(Cursor::new(buf));
//         // let mut parsing_template = ParsingTemplate::new();
//         // parsing_template.set_all();
//         // let mut reader = Reader::new(in_cursor, parsing_template).unwrap();
//         // let mut records = reader.records();
//         // let mut it = raw_records.iter();
//         // while let Some(rec) = records.next_rec() {
//         //     let rec_orig = it.next().unwrap();
//         //     let orig_map_q = rec_orig.get_bytes(&Fields::Mapq)[0];
//         //     let orig_pos = rec_orig
//         //         .get_bytes(&Fields::Pos)
//         //         .read_i32::<LittleEndian>()
//         //         .unwrap();
//         //     assert_eq!(rec.pos.unwrap(), orig_pos);
//         //     assert_eq!(rec.mapq.unwrap(), orig_map_q);
//         // }
//     }
// }
