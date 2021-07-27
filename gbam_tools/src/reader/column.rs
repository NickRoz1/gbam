use std::{
    borrow::Borrow,
    cell::RefCell,
    io::{Read, Result, Seek},
    ops::{Deref, DerefMut},
    sync::{Arc, Mutex},
};

use super::reader::{GbamRecord, ReadSeekSend};

use bam_tools::record::fields::Fields;
use lz4::block;
use memmap::Mmap;

use crate::meta::FileMeta;

// MAKE A TRAIT COLUMN AND THEN TWO SPECIALZIATIONS FOR FIXED SIZE FIELDS AND
// VARIABLE SIZED FIELDS (THE LATTER ONE REQUIRES OTHER INDEX FIELDS TO BE
// LOADED SIMULTANEOUSLY).

// Contains fields needed both for fixed sized fields and variable sized fields.
struct Inner<'a> {
    meta: &'a FileMeta,
    range_begin: usize,
    range_end: usize,
    field: Fields,
    block_num: usize,
    buffer: Vec<u8>,
    reader: RefCell<Mmap>,
}
/// Defines how columns will operate. It is needed since variable sized fields
/// columns also require parsing of additional fixed sized fields columns.
trait Column {
    // Fills GbamRecord field with data from corresponding BAM record.
    fn fill_record_field(&mut self, item_num: usize, rec: &mut GbamRecord) -> ();
}
/// GBAM file column. Responsible for fetching data.
struct FixedColumn<'a>(Inner<'a>);

impl<'a> Deref for FixedColumn<'a> {
    type Target = Inner<'a>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<'a> DerefMut for FixedColumn<'a> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl<'a> Column for FixedColumn<'a> {
    fn fill_record_field(&mut self, item_num: usize, rec: &mut GbamRecord) -> () {
        if let Some(block_num) = self.find_block(item_num) {
            fetch_block(&mut self.0, block_num);
            let block_len = self.meta.view_blocks_sizes(&self.field)[0] as usize;
            self.range_begin = block_num * block_len;
            self.range_end = self.range_begin + block_len;
        }
        rec.parse_from_bytes(&self.field, self.get_item(item_num));
    }
}

impl<'a> FixedColumn<'a> {
    // Finds blocks where record is located. None is returned if block is already loaded.
    fn find_block(&self, item_num: usize) -> Option<usize> {
        if item_num >= self.range_begin && item_num < self.range_end {
            return None;
        }
        // All blocks sizes are equal except maybe last one (since it's a fixed sized column).
        let block_len = self.meta.view_blocks_sizes(&self.field)[0];
        Some(item_num / block_len as usize)
    }

    fn get_item(&self, item_num: usize) -> &[u8] {
        let item_size = self.meta.get_field_size(&self.field).unwrap() as usize;
        let offset = item_num * item_size;
        &self.buffer[offset..offset + item_size]
    }
}

/// Column managing access to variable sized data. Utilized another column (for fixed sized fields) to index data.
struct VariableColumn<'a>(Inner<'a>, FixedColumn<'a>);

impl<'a> Column for VariableColumn<'a> {
    fn fill_record_field(&mut self, item_num: usize, rec: &mut GbamRecord) -> () {
        if let Some(block_num) = self.find_block(item_num) {
            fetch_block(block_num);
            let block_len = self.meta.view_blocks_sizes(&self.field)[0] as usize;
            self.range_begin = block_num * block_len;
            self.range_end = self.range_begin + block_len;
        }
        rec.parse_from_bytes(&self.field, self.get_item(item_num));
    }
}

fn fetch_block(inner_column: &mut Inner, block_num: usize) -> Result<()> {
    let field = &inner_column.field;
    let block_meta = inner_column.meta.view_blocks(field).get(block_num).unwrap();
    let reader = inner_column.reader.borrow_mut();
    let block_size = inner_column.meta.view_blocks_sizes(field)[block_num as usize];

    inner_column.buffer.clear();
    let mut data =
        &reader[(block_meta.seekpos as usize)..(block_meta.seekpos + block_size as u64) as usize];

    data.read_to_end(&mut inner_column.buffer)?;
    inner_column.block_num = block_num;

    Ok(())
}
