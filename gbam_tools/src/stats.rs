// pub type StatsComparator = fn(&[u8], &[u8]) -> Ordering;

// pub struct StatsCollector {
//     pub max_value: Vec<Stat>,
//     pub min_value: Vec<Stat>,
//     comparator: StatsComparator,
// }

// impl StatsCollector {
//     pub fn new(field: Fields, comparator: StatsComparator) -> Self {
//         Self {
//             field,
//             max_value: None,
//             min_value: None,
//             comparator,
//         }
//     }

//     pub fn update(&mut self, value: &[u8]) {
//         if self.max_value.is_none()
//             || (self.comparator)(value, &self.max_value.as_ref().unwrap()) == Ordering::Greater
//         {
//             self.max_value = Some(value.to_vec());
//         }
//         if self.min_value.is_none()
//             || (self.comparator)(value, &self.min_value.as_ref().unwrap()) == Ordering::Less
//         {
//             self.min_value = Some(value.to_vec());
//         }
//     }

//     pub fn reset(&mut self) {
//         self.max_value = None;
//         self.min_value = None;
//     }
// }

// pub fn refid_comparator(mut left: &[u8], mut right: &[u8]) -> Ordering {
//     let left_val = left.read_i32::<LittleEndian>().unwrap();
//     let right_val = right.read_i32::<LittleEndian>().unwrap();

//     left_val.cmp(&right_val)
// }
