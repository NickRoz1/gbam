use std::{
    borrow::Borrow,
    cell::RefCell,
    collections::{BTreeMap, BTreeSet},
    io::Result,
    rc::Rc,
    sync::Arc,
};

use super::record::GbamRecord;

use bam_tools::record::fields::Fields;
use byteorder::{LittleEndian, ReadBytesExt};
use flate2::write::GzDecoder;
use lz4::Decoder;
use memmap2::Mmap;

use crate::{meta::FileMeta, Codecs};

// Contains fields needed both for fixed sized fields and variable sized fields.
pub struct Inner {
    /// Arc is needed since this struct should work with PyO3 which sends struct between threads.
    meta: Arc<FileMeta>,
    range_begin: usize,
    range_end: usize,
    field: Fields,
    buffer: Vec<u8>,
    reader: Arc<Mmap>,
}

impl Inner {
    pub(crate) fn new(meta: Arc<FileMeta>, field: Fields, reader: Arc<Mmap>) -> Self {
        println!(
            "Field: {:?} block_num: {}",
            field,
            meta.view_blocks_sizes(&field)[0]
        );
        Inner {
            meta,
            range_begin: 0,
            range_end: 0,
            field,
            buffer: Vec::<u8>::new(),
            reader,
        }
    }
}

/// Defines how columns will operate. It is needed since variable sized fields
/// columns also require parsing of additional fixed sized fields columns.
pub trait Column {
    // Fills GbamRecord field with data from corresponding BAM record.
    fn fill_record_field(&mut self, item_num: usize, rec: &mut GbamRecord) -> ();
}

/// GBAM file column. Responsible for fetching data.
pub struct FixedColumn(Inner);

impl Column for FixedColumn {
    /// Fetches data into provider record buffer. If item is located outside of
    /// currently loaded data block, the new block will be loaded and
    /// decompressed.
    fn fill_record_field(&mut self, item_num: usize, rec: &mut GbamRecord) -> () {
        rec.parse_from_bytes(&self.0.field.clone(), self.get_item(item_num));
    }
}

impl FixedColumn {
    pub fn new(inner: Inner) -> Self {
        Self(inner)
    }
    fn get_item(&mut self, item_num: usize) -> &[u8] {
        if let Some(block_num) = self.find_block(item_num) {
            update_buffer(&mut self.0, block_num);
        }
        let item_size = self.0.meta.get_field_size(&self.0.field).unwrap() as usize;
        let offset = item_num * item_size;
        &self.0.buffer[offset..offset + item_size]
    }
    // Finds blocks where record is located. None is returned if block is already loaded.
    fn find_block(&self, item_num: usize) -> Option<usize> {
        if item_num >= self.0.range_begin && item_num < self.0.range_end {
            return None;
        }
        // All blocks sizes are equal except maybe last one (since it's a fixed sized column).
        let block_len = self.0.meta.view_blocks(&self.0.field)[0].numitems;
        // println!(
        //     "Field {:?} Item_num {}, divisor {}, result {}, len {}",
        //     self.0.field,
        //     item_num,
        //     block_len,
        //     item_num / block_len as usize,
        //     self.0.meta.view_blocks_sizes(&self.0.field).len()
        // );
        Some(item_num / block_len as usize)
    }
}

/// Column managing access to variable sized data. Utilizes another column (for fixed sized fields) to index data.
pub struct VariableColumn {
    inner: Inner,
    index: FixedColumn,
    // Used to quickly determine what block record belongs to.
    blocks: BTreeMap<usize, usize>,
}

impl Column for VariableColumn {
    fn fill_record_field(&mut self, item_num: usize, rec: &mut GbamRecord) -> () {
        rec.parse_from_bytes(&self.inner.field.clone(), self.get_item(item_num));
    }
}

impl VariableColumn {
    pub fn new(inner: Inner, index: FixedColumn) -> Self {
        Self {
            blocks: inner
                .meta
                .view_blocks(&inner.field)
                .iter()
                .enumerate()
                // Prefix sum.
                .scan(0, |acc, (count, x)| {
                    *acc = *acc + x.numitems;
                    Some((*acc as usize, count))
                })
                .collect(),
            inner,
            index,
        }
    }

    fn get_item(&mut self, item_num: usize) -> &[u8] {
        if let Some(block_num) = self.find_block(item_num) {
            update_buffer(&mut self.inner, block_num);
        }
        let mut read_offset =
            |n| self.index.get_item(n).read_u32::<LittleEndian>().unwrap() as usize;
        let start = match item_num {
            0 => 0,
            _ => read_offset(item_num - 1),
        };
        let end = read_offset(item_num);
        &self.inner.buffer[start..end]
    }

    // Finds blocks where record is located. None is returned if block is already loaded.
    fn find_block(&self, item_num: usize) -> Option<usize> {
        if item_num >= self.inner.range_begin && item_num < self.inner.range_end {
            return None;
        }
        // To determine what block record N is in, count elements smaller than N in the BTree.
        Some(
            self.blocks
                .range(..item_num)
                .next_back()
                .map_or(0, |(_, &val)| val),
        )
    }
}

fn update_buffer(inner: &mut Inner, block_num: usize) {
    println!("Doing heavy work");
    fetch_block(inner, block_num).unwrap();
    let block_len = inner.meta.view_blocks(&inner.field)[block_num].numitems;
    /// FIX HERE
    inner.range_begin = block_num * block_len;
    inner.range_end = inner.range_begin + block_len;
}

/// Fetch and decompress a data block.
fn fetch_block(inner_column: &mut Inner, block_num: usize) -> Result<()> {
    let field = &inner_column.field;
    let block_meta = inner_column.meta.view_blocks(field).get(block_num).unwrap();
    let reader = &inner_column.reader;
    let block_size = inner_column.meta.view_blocks_sizes(field)[block_num as usize];

    let data =
        &reader[(block_meta.seekpos as usize)..(block_meta.seekpos + block_size as u64) as usize];

    inner_column.buffer.clear();
    let codec = inner_column.meta.get_field_codec(field);

    decompress_block(data, &mut inner_column.buffer, codec).expect("Decompression failed.");
    Ok(())
}

fn decompress_block(source: &[u8], dest: &mut Vec<u8>, codec: &Codecs) -> std::io::Result<()> {
    use std::io::Write;
    match codec {
        Codecs::Gzip => {
            let mut decoder = GzDecoder::new(dest);
            decoder.write_all(source).unwrap();
            decoder.try_finish().unwrap();
        }
        Codecs::Lz4 => {
            let mut decoder = Decoder::new(source)?;
            std::io::copy(&mut decoder, dest).unwrap();
        }
    };
    Ok(())
}
