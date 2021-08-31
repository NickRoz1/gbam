use std::cmp::Ordering;

use bam_tools::record::{bamrawrecord::BAMRawRecord, fields::Fields};

use byteorder::{LittleEndian, ReadBytesExt, WriteBytesExt};

struct StatsCollector<T>
where
    T: Clone,
{
    field: Fields,
    max_value: Option<T>,
    min_value: Option<T>,
    comparator: Box<dyn Fn(&T, &T) -> Ordering>,
}

impl<T> StatsCollector<T>
where
    T: Clone,
{
    pub fn new(max_value: Option<T>, min_value: Option<T>, field: Fields) -> Self {
        Self {
            field,
            max_value,
            min_value,
            comparator: get_comparator(&field),
        }
    }

    pub fn update(&mut self, value: &T) {
        if self.max_value.is_none()
            || (self.comparator)(value, &self.max_value.as_ref().unwrap()) == Ordering::Greater
        {
            self.max_value = Some(value.clone());
        }
        if self.min_value.is_none()
            || (self.comparator)(value, &self.min_value.as_ref().unwrap()) == Ordering::Less
        {
            self.min_value = Some(value.clone());
        }
    }
}

pub(crate) trait CollectStats<T> {
    fn collect_stats(&mut self, rec: &mut BAMRawRecord);
}

impl<T> CollectStats<T> for StatsCollector<T>
where
    T: Clone,
    StatsCollector<T>: ConvertFromBytes<T>,
{
    fn collect_stats(&mut self, rec: &mut BAMRawRecord) {
        self.update_with_rec(&rec);
    }
}

impl<T> StatsCollector<T>
where
    T: Clone,
    StatsCollector<T>: ConvertFromBytes<T>,
{
    pub fn update_with_rec(&mut self, rec: &BAMRawRecord) {
        self.update(&Self::from_bytes(rec.get_bytes(&self.field)))
    }
}

fn get_comparator<T>(field: &Fields) -> Box<dyn Fn(&T, &T) -> Ordering> {
    Box::new(|a, b| Ordering::Equal)
}

impl<T: Clone> StatsCollector<T>
where
    StatsCollector<T>: ConvertToBytes<T>,
{
    fn get_max_bytes(&self) -> Option<Vec<u8>> {
        self.max_value
            .as_ref()
            .and_then(|data| Some(Self::into_bytes(data.into())))
    }
    fn get_min_bytes(&self) -> Option<Vec<u8>> {
        self.min_value
            .as_ref()
            .and_then(|data| Some(Self::into_bytes(data.into())))
    }
}

pub(crate) trait ConvertToBytes<T> {
    fn into_bytes(val: &T) -> Vec<u8>;
}

impl ConvertToBytes<i32> for StatsCollector<i32> {
    fn into_bytes(val: &i32) -> Vec<u8> {
        let mut wtr = vec![];
        wtr.write_i32::<LittleEndian>(*val).unwrap();
        wtr
    }
}

impl ConvertToBytes<u16> for StatsCollector<u16> {
    fn into_bytes(val: &u16) -> Vec<u8> {
        let mut wtr = vec![];
        wtr.write_u16::<LittleEndian>(*val).unwrap();
        wtr
    }
}

impl ConvertToBytes<u8> for StatsCollector<u8> {
    fn into_bytes(val: &u8) -> Vec<u8> {
        vec![*val]
    }
}

pub(crate) trait ConvertFromBytes<T> {
    fn from_bytes(val: &[u8]) -> T;
}

impl ConvertFromBytes<i32> for StatsCollector<i32> {
    fn from_bytes(data: &[u8]) -> i32 {
        (&data[..]).read_i32::<LittleEndian>().unwrap()
    }
}

impl ConvertFromBytes<u16> for StatsCollector<u16> {
    fn from_bytes(data: &[u8]) -> u16 {
        (&data[..]).read_u16::<LittleEndian>().unwrap()
    }
}

impl ConvertFromBytes<u8> for StatsCollector<u8> {
    fn from_bytes(data: &[u8]) -> u8 {
        data[0]
    }
}
