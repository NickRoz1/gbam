use std::cmp::Ordering;

use bam_tools::record::fields::Fields;

pub struct StatsCollector {
    field: Fields,
    pub max_value: Option<Vec<u8>>,
    pub min_value: Option<Vec<u8>>,
    comparator: Box<dyn Fn(&[u8], &[u8]) -> Ordering>,
}

impl StatsCollector {
    pub fn new(field: Fields, comparator: Box<dyn Fn(&[u8], &[u8]) -> Ordering>) -> Self {
        Self {
            field,
            max_value: None,
            min_value: None,
            comparator,
        }
    }

    pub fn update(&mut self, value: &[u8]) {
        if self.max_value.is_none()
            || (self.comparator)(value, &self.max_value.as_ref().unwrap()) == Ordering::Greater
        {
            self.max_value = Some(value.to_vec());
        }
        if self.min_value.is_none()
            || (self.comparator)(value, &self.min_value.as_ref().unwrap()) == Ordering::Less
        {
            self.min_value = Some(value.to_vec());
        }
    }

    pub fn reset(&mut self) {
        self.max_value = None;
        self.min_value = None;
    }
}
