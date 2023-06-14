use bam_tools::MEGA_BYTE_SIZE;
use bam_tools::record::fields::Fields;
use byteorder::ReadBytesExt;
use std::arch::asm;
use std::cmp::min;
use std::convert::TryInto;
use std::io::{Write, BufWriter, StdoutLock, Seek, Read, BufReader};
use std::ops::{RangeInclusive, Range};
use std::sync::Arc;
use std::{cmp::Ordering, collections::HashMap, time::Instant};
use std::fs::File;
use crate::query::cigar::Cigar;
use crate::reader::parse_tmplt::ParsingTemplate;
use crate::reader::record::parse_cigar;
use crate::utils::bed;
/// This module provides function for fast querying of read depth.
use crate::meta::{BlockMeta, FileMeta};
use crate::reader::{reader::Reader, record::GbamRecord};
use std::path::{PathBuf};
use crossbeam::channel::{Receiver, Sender, bounded};
use std::thread;
use std::thread::JoinHandle;
use super::int2str::{i32toa_countlut, u32toa_countlut};
use flate2::Compression;
use flate2::write::GzEncoder;
use std::io::SeekFrom;
use crate::reader::column::decompress_block;

#[allow(dead_code)]
fn panic_err() {
    panic!("The query you entered is incorrect. The format is as following: <ref name>:<position>\ne.g. chr1:1257\n");
}

fn process_range(preparsed_records: Arc<Vec<DepthUnit>>, index_file: Option<Arc<Vec<u32>>>, rec_range: Range<usize>, mut scan_line: Vec<i32>, target_id: i32) -> Vec<i32> {
    // let mut rec = GbamRecord::default();
    for idx in rec_range {
        let rec = preparsed_records[index_file.as_ref().unwrap()[idx] as usize];
        if rec.refid != target_id {
            break;
        }
        if rec.cigar == 0 {
            continue;
        }
        let read_start: usize = rec.pos as usize;
        let mut base_cov = rec.cigar as usize;
        if(rec.flag & 0b11100000100) != 0 {
            base_cov = 0;
        }
        if base_cov == 0 {
            base_cov = 1;
        }
      
        let read_end = read_start + base_cov;

        scan_line[read_start] += 1;
        scan_line[read_end] -= 1;
        // buf.increments.push(read_start);
        // buf.decrements.push(read_end);
    }
    scan_line
}

fn calc_depth(preparsed_records: Arc<Vec<DepthUnit>>, file_meta: Arc<FileMeta>, index_file: Option<Arc<Vec<u32>>>, number_of_records: usize, ref_id: i32, mut coverage_arr: Vec<i32>, ref_len: usize) -> Vec<i32> {
    coverage_arr.resize(ref_len+1, 0);

    // let lower_bound = if let Some(block_num) = find_leftmost_block(ref_id, file_meta.view_blocks(&Fields::RefID)) {
    //     block_num as usize
    // }
    // else {
    //     // This refid was not found in any block
    //     return coverage_arr;
    // };
    // let upper_bound = find_rightmost_block(ref_id, file_meta.view_blocks(&Fields::RefID)) as usize;
    let mut first_rec:i64 = -1;
    
    let mut last_rec:i64=  preparsed_records.len() as i64;
    
  
    while(last_rec - first_rec > 1){
        let mid: usize = ((first_rec + last_rec)/2) as usize;
        let buf = preparsed_records[index_file.as_ref().unwrap()[mid] as usize];
        if buf.refid >= ref_id || buf.refid == -1 {
            last_rec = mid as i64;
        }
        else{
            first_rec = mid as i64;
        }
    }
    first_rec += 1;
    let amount = preparsed_records.len();

    if first_rec as usize == amount || preparsed_records[index_file.as_ref().unwrap()[first_rec as usize] as usize].refid != ref_id {
        return coverage_arr;
    }


    // let last_rec = std::cmp::min(upper_bound*file_meta.view_blocks(&Fields::RefID)[0].numitems as usize, number_of_records-1);

    // let mut temp_reader = Reader::new_with_meta(gbam_file.try_clone().unwrap(), ParsingTemplate::new_with(&[Fields::RefID, Fields::Pos, Fields::RawCigar]), &file_meta).unwrap();
    // let mut rec = GbamRecord::default();
    // temp_reader.fill_record(last_rec, &mut rec);

    // let read_start: usize = rec.pos.unwrap().try_into().unwrap();
    // let base_cov = rec.cigar.as_ref().unwrap().base_coverage() as usize;
    // let read_end = read_start + base_cov;

    // Loads of page faults here.
    

    // dbg!("Allocated {}", ref_len);

    let mut coverage = process_range(preparsed_records, index_file, first_rec as usize..amount, coverage_arr, ref_id);
    let mut acc = 0;
    for slot in coverage.iter_mut() {
        acc += *slot;
        *slot = acc; 
    }
    coverage
}

#[derive(Default, Clone, Copy)]
struct DepthUnit {
    refid: i32,
    pos: i32,
    cigar: u32,
    flag: u16,
}

pub fn main_depth(mut gbam_file: File, bed_file: Option<&PathBuf>, index_file: Option<Arc<Vec<u32>>>, bed_cli_request: Option<String>, _mapq: Option<u32>, bed_gz_path: Option<PathBuf>, thread_num: Option<usize>){    
    let mut queries = HashMap::<String, Vec<(u32, u32)>>::new();
    if let Some(bed_path) = bed_file {
        queries = bed::parse_bed_from_file(bed_path).expect("BED file is corrupted.");
    } 
    if let Some(query) = bed_cli_request {
        queries.extend(bed::parse_bed(&mut query.as_bytes()).unwrap().into_iter());
    }

    let mut reader = Reader::new(gbam_file.try_clone().unwrap(), ParsingTemplate::new()).unwrap();
    let file_meta = reader.file_meta.clone();
    let ref_seqs = file_meta.get_ref_seqs().clone();
    let chr_to_ref_id = get_chr_name_mapping(ref_seqs.iter().map(|(chr, _)| chr), &mut reader);
    let number_of_records = reader.amount;
    drop(reader);

    // Calculate for whole file.
    if queries.is_empty() {
        ref_seqs.iter().for_each(|(chr, len)| {queries.insert(chr.clone(), vec![(0, *len)]);});
    }

    let mut buffers = vec![Vec::<i32>::new()];
    if thread_num.is_some(){
        buffers = vec![Vec::<i32>::new();std::cmp::min(thread_num.unwrap(), 8)];
    }

    type VectorOfSendersAndReceivers= Vec::<Option<(Sender<(Arc<Vec<DepthUnit>>, Arc<FileMeta>, Option<Arc<Vec<u32>>>, usize, i32, Vec<i32>, usize, String)>,Receiver<(String, Vec<i32>)>)>>;
    let mut circular_buf_channels = VectorOfSendersAndReceivers::new();
    (0..buffers.len()).for_each(|_|circular_buf_channels.push(None));
    let mut handles: Vec::<JoinHandle<()>> = Vec::new();

    let mut idx = 0;
    // let mut coverage_arr: Vec<i64> = Vec::new(); 
    // coverage_arr.reserve(longest_chr as usize);

    
    let mut iter = ref_seqs.iter();
    let mut accum = 0;  
    let mut bed_gz_printer = bed_gz_path.map(BedGzPrinter::new);

    let st = std::io::stdout();
    let lock = st.lock();
    let mut printer = ConsolePrinter::new(lock);

    let mut temp: Reader = Reader::new_with_meta(gbam_file.try_clone().unwrap(), ParsingTemplate::new_with(&[Fields::RefID, Fields::Pos, Fields::Flags]), &file_meta, None).unwrap();
    dbg!(temp.amount);
    let mut preparsed = vec![DepthUnit::default();temp.amount];

    let mut it = temp.records();
    let mut var_i = 0;
    let mut now = Instant::now();
    
    while let Some(rec) = it.next_rec() {
     
        // if var_i % 100_000 == 0{
        //     dbg!("Processed records:");
        //     dbg!(var_i);
        // }
        preparsed[var_i].refid = rec.refid.unwrap();
        preparsed[var_i].pos = rec.pos.unwrap();
        // preparsed[var_i].cigar = rec.cigar.as_ref().unwrap().base_coverage();
        preparsed[var_i].flag = rec.flag.unwrap();
        var_i += 1;
    }   

    dbg!("Took to fill with fixed sized fields:");
    dbg!(now.elapsed().as_millis());
    
    now = Instant::now();
    
    
    let mut buf = Vec::new();
    let mut decompressed_buf = Vec::new();

    let codec = file_meta.get_field_codec(&Fields::NCigar);
    let mut indices = Vec::new();
    let mut read_manual = BufReader::with_capacity(MEGA_BYTE_SIZE, gbam_file.try_clone().unwrap());
    for block in file_meta.view_blocks(&Fields::NCigar){
        let available_in_block = block.numitems;
        buf.resize(block.block_size as usize, 0);
        decompressed_buf.resize(block.uncompressed_size as usize, 0);
        read_manual.seek(SeekFrom::Start(block.seekpos)).unwrap();
        read_manual.read_exact(&mut buf).unwrap();
        decompress_block(&buf, &mut decompressed_buf, codec).unwrap();
        let mut slice = &decompressed_buf[..];
        for _ in 0..available_in_block {
            // dbg!(slice.len());
            indices.push(slice.read_u32::<byteorder::LittleEndian>().unwrap() as usize);
        }
    }

    assert!(temp.amount == indices.len());
    let codec = file_meta.get_field_codec(&Fields::RawCigar);

    let mut read_in_total = 0;
    let mut cigar_buf = Cigar::new(Vec::new());
    var_i = 0;
    for block in file_meta.view_blocks(&Fields::RawCigar){
        let available_in_block = block.numitems;
        let mut read_currently = 0;
        buf.resize(block.block_size as usize, 0);
        decompressed_buf.resize(block.uncompressed_size as usize, 0);
        read_manual.seek(SeekFrom::Start(block.seekpos)).unwrap();
        read_manual.read_exact(&mut buf).unwrap();
        // decompress_block(&buf, &mut decompressed_buf, codec).unwrap();
        
        while read_currently < available_in_block {
            let the_cigar_end = indices[read_in_total];
            let mut the_cigar_begin = 0;
            if read_currently > 0 {
                the_cigar_begin = indices[read_in_total-1];
            }
            
            if var_i % 100_000 == 0{
                dbg!("Processed records:");
                dbg!(var_i);
            }
            // let actual_cigar = &decompressed_buf[the_cigar_begin..the_cigar_end];
            // parse_cigar(actual_cigar, &mut cigar_buf);

            preparsed[var_i].cigar = cigar_buf.base_coverage();
            var_i += 1;
            read_currently += 1;
            read_in_total += 1;
        }
        
    }

    let arc_of_records = Arc::new(preparsed);

    dbg!("Finished parsing all records to RAM buffer.");

    loop {
        // dbg!(buffers.len()); 
        if idx == circular_buf_channels.len() {
            idx = 0;
        }
        if circular_buf_channels[idx].is_some() {
            let (thread_chr, mut coverage_arr) = circular_buf_channels[idx].as_mut().unwrap().1.recv().unwrap();

            if let Some(bed_regions) = queries.get(&thread_chr) {
                // coverage_arr.resize(*ref_len as usize, 0);
                // let ref_id = chr_to_ref_id.get(chr).unwrap().unwrap();
                // buffers = calc_depth(gbam_file.try_clone().unwrap(), file_meta.clone(), number_of_records, ref_id, &mut coverage_arr, buffers);
    
                // 41641770854
                
                // printer.set_chr(thread_chr.clone());
                let now = Instant::now();
                if bed_gz_printer.is_none() {
                    
                    for bed_region in bed_regions {
                        let st = bed_region.0 as usize;
                        let en = min((bed_region.1) as usize, coverage_arr.len());
                        for coord in st..en {
                            unsafe {
                                if coverage_arr.get_unchecked(coord) > &0 {
                                    printer.write_efficient(&thread_chr, (coord) as u32, *coverage_arr.get_unchecked(coord));
                                }
                            }
                        }
                    }
                    
                }
                else {
                    
                    for bed_region in bed_regions {
                        let st = bed_region.0;
                        let en = min(bed_region.1, (coverage_arr.len()) as u32);
                        let mut prev_coord = None;
                        let mut prev_depth = None;
                        
                        for coord in st..en {
                            unsafe {
                                let cur_depth = *coverage_arr.get_unchecked(coord as usize);
                                // if cur_depth > 0 {
                                    if prev_coord.is_none(){
                                        prev_coord = Some(coord);
                                        prev_depth = Some(cur_depth);
                                    } 
                                    else if prev_depth.unwrap() != cur_depth {
                                        
                                            bed_gz_printer.as_mut().unwrap().write_region(&thread_chr, prev_coord.unwrap(), coord, prev_depth.unwrap());
                                            prev_depth = Some(cur_depth);
                                            prev_coord = Some(coord);
                                        
                                    }
                                // }
                            }
                        }

                        bed_gz_printer.as_mut().unwrap().write_region(&thread_chr, prev_coord.unwrap() , en , prev_depth.unwrap());                        
                    }
                }
                accum += now.elapsed().as_millis();
                coverage_arr.clear();
            }
            
            buffers.push(coverage_arr);
        }

        let next_chr = iter.next(); 
        
        if let Some((chr, ref_len)) = next_chr {
            if circular_buf_channels[idx].is_none() {
                let (s, r) = bounded(1);
                let (ready_s, ready_r) = bounded(1);
                let handle = thread::spawn(move || {
                    for task  in r {
                        let (preparsed, meta, index_file, number_of_records, ref_id, buf, t_ref_len, t_chr) = task;
                        ready_s.send((t_chr, calc_depth(preparsed, meta, index_file, number_of_records, ref_id, buf, t_ref_len))).unwrap();
                    } 
                });
                circular_buf_channels[idx] = Some((s, ready_r));
                handles.push(handle);
            }

            let ref_id = chr_to_ref_id.get(chr).unwrap().unwrap();
            let buf = buffers.pop().unwrap();
            let meta = file_meta.clone();
            let file = gbam_file.try_clone().unwrap();
            let index = index_file.as_ref().map(|f| f.clone());
            let t_chr = chr.clone();
            let t_ref_len = *ref_len as usize;
            circular_buf_channels[idx].as_mut().unwrap().0.send((arc_of_records.clone(), meta, index, number_of_records, ref_id, buf, t_ref_len, t_chr)).unwrap();
        }

        if buffers.len() == circular_buf_channels.len(){
            break;
        }

        idx += 1;
    }

    // Drop all senders so that threads get closed.
    circular_buf_channels.clear();

    for h in handles {
        h.join().unwrap();
    }

    dbg!(accum);
    // Shouldn't allocate more.
    // assert!(coverage_arr.capacity() == longest_chr as usize);
}

fn get_chr_name_mapping<'a, I>(ref_ids: I, reader: &mut Reader) -> HashMap<String, Option<i32>>
where
    I: Iterator<Item = &'a String>,
{
    let mut name_to_ref_id = HashMap::<String, i32>::new();
    reader
        .file_meta
        .get_ref_seqs()
        .iter()
        .enumerate()
        .for_each(|(i, (name, _))| {
            name_to_ref_id.insert(name.clone(), i as i32);
        });
    // https://github.com/biod/sambamba/blob/3eff9a2d8bb3097b92c72752be3c6b42dd1c59b7/BioD/bio/std/hts/bam/read.d#L902
    ref_ids
        .map(|id| (id.to_owned(), name_to_ref_id.get(id).cloned()))
        .collect::<HashMap<String, Option<i32>>>()
}


// chrM    15268   438
// chrM    15269   439
// chrM    15270   425
// chrM    15271   420
// chrM    15272   418
// chrM    15273   407
// chrM    15274   268
// chrM    15275   278
// chrM    15276   281
// Approach as in https://github.com/brentp/mosdepth

fn find_leftmost_block(id: i32, block_metas: &Vec<BlockMeta>) -> Option<i64> {
    let mut left: i64 = -1;
    let mut right: i64 = block_metas.len() as i64;
    while (right-left) > 1 {
        let mid = (left + right) / 2;
        let max_val = &block_metas[mid as usize].stats.as_ref().unwrap().max_value;
        match max_val.cmp(&id) {
            Ordering::Equal | Ordering::Greater => right = mid,
            Ordering::Less => left = mid,
        }
    }
    if right as usize == block_metas.len() || block_metas[right as usize].stats.as_ref().unwrap().min_value > id {
        return None;
    }
    Some(right)
}

fn find_rightmost_block(id: i32, block_metas: &Vec<BlockMeta>) -> i64 {
    let mut left: i64 = -1;
    let mut right: i64 = block_metas.len() as i64;
    while (right-left) > 1 {
        let mid = (left + right) / 2;
        let min_val = &block_metas[mid as usize].stats.as_ref().unwrap().min_value;
        match min_val.cmp(&id) {
            Ordering::Equal | Ordering::Less => left = mid,
            Ordering::Greater => right = mid,
        }
    }
    right
}


struct ConsolePrinter<'a> {
    buffer: [u8; 400],
    stdout: BufWriter<StdoutLock<'a>>
 
}
impl<'a> ConsolePrinter<'a> {
    pub fn new(stdout_lock: StdoutLock<'a>) -> Self {
        let stdout = BufWriter::with_capacity(64 * 1024, stdout_lock);
        Self {  
            buffer: [0;400],
            stdout,
        }
    }

    /// Done in reversed direction because we don't know what is the size of integers beforehand.
    pub fn write_efficient(&mut self, reversed_chr: &str,  coord: u32,  depth: i32){
        let mut buff_ptr = self.buffer.as_mut_ptr();
        let orig: *mut u8 = self.buffer.as_mut_ptr();
        unsafe {
            for &ch in reversed_chr.as_bytes() {
                *buff_ptr = ch;
                buff_ptr = buff_ptr.add(1);
            }
            *buff_ptr = b'\t';
            buff_ptr = buff_ptr.add(1);
            buff_ptr = u32toa_countlut(coord, buff_ptr);
            *buff_ptr = b'\t';
            buff_ptr = buff_ptr.add(1);
            buff_ptr = i32toa_countlut(depth, buff_ptr);
            *buff_ptr = b'\n';
            buff_ptr = buff_ptr.add(1);
            self.stdout.write_all(&self.buffer[..(buff_ptr as usize - orig as usize)]).unwrap();
        }
    }
}

// chr1    0       10571   0
// chr1    10571   10598   1
// chr1    10598   15904   0
// chr1    15904   15931   1
// chr1    15931   16375   0
// chr1    16375   16402   1
// chr1    16402   16579   0
// chr1    16579   16606   4
// chr1    16606   18816   0
// chr1    18816   18843   1
// chr1    18843   19754   0
// chr1    19754   19781   1
struct BedGzPrinter{
    buffer: [u8; 400],
    compressor: GzEncoder<BufWriter<File>>,

}
impl BedGzPrinter {
    pub fn new(path: PathBuf) -> Self {
        let  file = File::create(path).expect("Failed to create depth file.");
        Self {  
            buffer: [0;400],
            compressor: GzEncoder::new(BufWriter::with_capacity(64 * 1024, file), Compression::default()),
        }
    }

    /// Done in reversed direction because we don't know what is the size of integers beforehand.
    pub fn write_region(&mut self, chr: &str, prev_coord: u32, coord: u32, prev_depth: i32){
        let mut buff_ptr = self.buffer.as_mut_ptr();
        let orig: *mut u8 = self.buffer.as_mut_ptr();
        unsafe {
            for &ch in chr.as_bytes() {
                *buff_ptr = ch;
                buff_ptr = buff_ptr.add(1);
            }
            *buff_ptr = b'\t';
            buff_ptr = buff_ptr.add(1);
            buff_ptr = u32toa_countlut(prev_coord, buff_ptr);
            *buff_ptr = b'\t';
            buff_ptr = buff_ptr.add(1);
            buff_ptr = u32toa_countlut(coord, buff_ptr);
            *buff_ptr = b'\t';
            buff_ptr = buff_ptr.add(1);
            buff_ptr = i32toa_countlut(prev_depth, buff_ptr);
            *buff_ptr = b'\n';
            buff_ptr = buff_ptr.add(1);
            self.compressor.write_all(&self.buffer[..(buff_ptr as usize - orig as usize)]).unwrap();
        }
    }
}