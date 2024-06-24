
use std::path::{Path, PathBuf};
use std::fs::File;
use std::io::{Write, Read};
use std::io::BufReader;
use std::cmp::Ordering;
use byteorder::{LittleEndian, ReadBytesExt, WriteBytesExt};
use std::collections::BinaryHeap;
use std::io::Cursor;

use crate::reader::record::GbamRecord;
extern crate rand;

use rand::Rng;
use rayon::slice::ParallelSliceMut;

const ONE_GB : usize = 1e9 as usize;
type Comparator = fn(&GbamRecord, &GbamRecord) -> Ordering;
type MergeRecord = (GbamRecord, File);

fn picard_compare(lhs: &GbamRecord, rhs: &GbamRecord) -> Ordering {

    let mut compareDifference  = 0;

    if compareDifference == 0 {
        compareDifference = lhs.refid.unwrap() - rhs.refid.unwrap();
    }
    if compareDifference == 0 {
        compareDifference = lhs.pos.unwrap() - rhs.pos.unwrap();
    }
    if compareDifference == 0 {
        compareDifference = lhs.is_reverse() as i32 - rhs.is_reverse() as i32;
    }
    if compareDifference == 0 {
        compareDifference = lhs.read2ReferenceIndex - rhs.read2ReferenceIndex;
    }
    if compareDifference == 0 {
        compareDifference = lhs.read2Coordinate - rhs.read2Coordinate;
    }

    if compareDifference == 0 {
        compareDifference = lhs.getTile() - rhs.getTile();
    }

    if compareDifference == 0 {
        compareDifference = lhs.getX() - rhs.getX();
    }

    if compareDifference == 0 {
        compareDifference = lhs.getY() - rhs.getY();
    }

    // The following is arbitrary and is only included for completeness.
    // Other implementations may chose to forgo this tiebreak if they do not have
    // access to the index-in-file of the records (e.g. SPARK implmentations)

    if compareDifference == 0 {
        compareDifference = Long.compare(lhs.read1IndexInFile, rhs.read1IndexInFile);
    }
    if compareDifference == 0 {
        compareDifference = Long.compare(lhs.read2IndexInFile, rhs.read2IndexInFile);
    }

}

struct Wrapper(GbamRecord, Comparator, BufReader<File>);

impl Ord for Wrapper {
    fn cmp(&self, other: &Self) -> Ordering {
        // Smallest go first.
        self.1(&other.0, &self.0)
    }
}


impl PartialOrd for Wrapper {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        panic!("WH");
    }
}

impl Eq for Wrapper {}

impl PartialEq for Wrapper {
    fn eq(&self, other: &Self) -> bool {
        panic!("WH");
    }
}

pub struct SortedStorageIter {
    sorted_storage: SortedStorage,
    heap: BinaryHeap<Wrapper>,
}


fn read_rec(cur: &mut BufReader<File>) -> Option<GbamRecord> {
    if let Ok(len) = cur.read_u64::<LittleEndian>() {
        let mut buf = Vec::<u8>::new();
        buf.resize(len as usize, 0);
        cur.read_exact(&mut buf).unwrap();
        bincode::deserialize(&buf).unwrap()
    }
    else {
        return None
    } 
}

impl SortedStorageIter {
    pub fn new(mut sorted_storage: SortedStorage) -> SortedStorageIter {
        let mut heap = BinaryHeap::new();
        for f in sorted_storage.temp_files.drain(..) {
            let mut reader = BufReader::new(f);
            heap.push(Wrapper(read_rec(&mut reader).unwrap(), sorted_storage.comp, reader));
        }
        SortedStorageIter {
            sorted_storage,
            heap
        }
    }
}
impl Iterator for SortedStorageIter {
    type Item = GbamRecord;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(mut wrap) = self.heap.pop(){
            let rec = wrap.0;
            self.heap.push(Wrapper(read_rec(&mut wrap.2).unwrap(), wrap.1, wrap.2));
            Some(rec)
        }
        else{
            None
        }

    }
}

pub struct SortedStorage {
    storage: Vec<GbamRecord>,
    temp_files: Vec<File>,
    comp: Comparator
}

// If there are more than 1GB of records, spills the rest into temp directory.
impl SortedStorage {
    pub fn sorted_storage(comp: Comparator) -> SortedStorage {
        SortedStorage {
            comp,
            temp_files: Vec::new(),
            storage: Vec::new(),
        }
    }

    pub fn add(&mut self, rec: GbamRecord){
        if self.storage.len() >= ONE_GB {
            self.flush();
        }
        self.storage.push(rec);
    }

    fn flush(&mut self){
        let dir =  std::env::temp_dir();
        let new_name =  md5::compute(rand::thread_rng().gen_range(33..(1e9 as u32)).to_string().as_bytes());
        let mut file = File::create(dir.as_path().join(format!("{:x}", new_name))).unwrap();

        self.storage[..].par_sort_by(self.comp);

        for r in self.storage.iter() {
            let d = bincode::serialize(&r).unwrap();
            file.write_u64::<LittleEndian>(d.len() as u64).unwrap();
            file.write_all(&d);
        }
    
        self.temp_files.push(file);
        self.storage.clear();
    }
}