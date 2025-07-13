use super::{reader::Reader, record::GbamRecord};

/// Iterates over GBAM file.
pub struct Records<'a> {
    reader: &'a mut Reader,
    cur_rec: usize,
    rec_amount: usize,
    buf: GbamRecord,
}

impl<'a> Records<'a> {
    pub fn new(reader: &'a mut Reader) -> Self {
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
