use super::{reader::Reader, record::GbamRecord};

/// Iterates over GBAM file.
pub struct Records<'a> {
    reader: &'a mut Reader,
    cur_rec: usize,
    rec_amount: usize,
    buf: GbamRecord,
}

impl<'a> Records<'a> {
    pub fn new(reader: &'a mut Reader, cur_rec: usize, rec_amount: usize, buf: GbamRecord) -> Self {
        Self {
            reader,
            cur_rec,
            rec_amount,
            buf,
        }
    }

    pub fn next_rec(&mut self) -> Option<&GbamRecord> {
        if self.cur_rec == self.rec_amount {
            return None;
        }
        for &field in self.reader.parsing_template.get_active_fields_iter() {
            self.reader.columns[field as usize]
                .as_mut()
                .unwrap()
                .fill_record_field(self.cur_rec, &mut self.buf);
        }
        Some(&self.buf)
    }
}
