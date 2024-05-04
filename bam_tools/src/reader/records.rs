// Source: https://github.com/zaeleus/noodles/blob/master/noodles-bam/src/reader/records.rs

use std::io::{self};

use super::Reader;

/// An iterator over records of a BAM reader.
///
/// This is created by calling [`Reader::records`].
pub struct Records<'a> {
    reader: &'a mut Reader,
    record: Vec<u8>,
}

impl<'a> Records<'a> {
    pub(crate) fn new(reader: &'a mut Reader) -> Records<'_> {
        Self {
            reader,
            record: Vec::default(),
        }
    }

    pub fn next_rec(&mut self) -> Option<io::Result<&Vec<u8>>> {
        match self.reader.read_record(&mut self.record) {
            Ok(0) => None,
            Ok(_) => Some(Ok(&self.record)),
            Err(e) => Some(Err(e)),
        }
    }
}

impl<'a> Iterator for Records<'a> {
    type Item = io::Result<Vec<u8>>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.reader.read_record(&mut self.record) {
            Ok(0) => None,
            Ok(_) => Some(Ok(self.record.clone())),
            Err(e) => Some(Err(e)),
        }
    }
}
