use super::{parse_tmplt::ParsingTemplate, reader::Reader, record::GbamRecord};
use std::fs::File;

#[cfg(feature = "python-ffi")]
use pyo3::prelude::*;

/// Iterates over GBAM file.
pub struct Records<'b> {
    reader: &'b mut Reader<'b>,
    cur_rec: usize,
    rec_amount: usize,
    buf: GbamRecord,
}

impl<'b> Records<'b> {
    pub fn new(reader: &'b mut Reader<'b>) -> Self {
        Self {
            rec_amount: reader.amount,
            reader,
            cur_rec: 0,
            buf: GbamRecord::default(),
        }
    }

    pub fn next_rec(&mut self) -> Option<&GbamRecord> {
        if self.cur_rec == self.rec_amount {
            return None;
        }
        self.reader.fill_record(self.cur_rec, &mut self.buf);
        self.cur_rec += 1;
        Some(&self.buf)
    }
}

// #[cfg_attr(feature = "python-ffi", pyclass)]
// /// This class was created because no generics (lifetimes) allowed in PyO3.
// pub struct PyRecords {
//     reader: Reader,
//     cur_rec: usize,
//     rec_amount: usize,
//     buf: GbamRecord,
// }

// impl PyRecords {
//     pub fn next_rec(&mut self) -> Option<&GbamRecord> {
//         if self.cur_rec == self.rec_amount {
//             return None;
//         }
//         self.reader.fill_record(self.cur_rec, &mut self.buf);
//         self.cur_rec += 1;
//         Some(&self.buf)
//     }
// }

// Methods to be called from Python.
#[cfg(feature = "python-ffi")]
#[pymethods]
impl PyRecords {
    #[allow(clippy::unnecessary_wraps)]
    fn next_record(&mut self) -> PyResult<Option<GbamRecord>> {
        let next = self.next_rec().cloned();
        Ok(next)
    }

    /// Create new reader for file at path
    #[new]
    pub fn new_reader(path: &str, tmplt: ParsingTemplate) -> Self {
        let file = File::open(path).unwrap();
        let reader = Reader::new(file, tmplt).unwrap();
        PyRecords {
            rec_amount: reader.amount,
            reader,
            cur_rec: 0,
            buf: GbamRecord::default(),
        }
    }
}
