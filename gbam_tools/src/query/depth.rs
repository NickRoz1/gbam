use bam_tools::record::fields::Fields;
use itertools::Itertools;
use std::ascii::AsciiExt;
use std::collections::BTreeMap;
use std::convert::TryInto;
use std::io::Write;
use std::ops::{Range, RangeInclusive};
use std::{cmp::Ordering, collections::HashMap, convert::TryFrom, time::Instant};

/// This module provides function for fast querying of read depth.
use crate::meta::{BlockMeta};
use crate::reader::{reader::generate_block_treemap, reader::Reader, record::GbamRecord};

type Region = RangeInclusive<u32>;

use byteorder::{LittleEndian, ReadBytesExt};

// TODO: Merge bed regions with tolerance to form super regions. Then do
// sweepline for each of the super regions. Print the depth results for each of
// the regions by taking chunks from the super regions or doing in process while
// calculating to avoid double iterating the super region sweep line array.

fn parse_bytes(bytes: &Vec<u8>) -> i32 {
    (&bytes[..]).read_i32::<LittleEndian>().unwrap()
}

fn get_chr_name_mapping<'a, I>(ref_ids: I, reader: &mut Reader) -> HashMap<String, Option<i32>>
where
    I: Iterator<Item = &'a String>,
{
    let mut name_to_ref_id = HashMap::<String, i32>::new();
    reader
        .file_meta
        .get_ref_seqs()
        .into_iter()
        .enumerate()
        .for_each(|(i, (name, _))| {
            name_to_ref_id.insert(name.clone(), i as i32);
        });
    // https://github.com/biod/sambamba/blob/3eff9a2d8bb3097b92c72752be3c6b42dd1c59b7/BioD/bio/std/hts/bam/read.d#L902
    ref_ids
        .map(|id| (id.to_owned(), name_to_ref_id.get(id).cloned()))
        .collect::<HashMap<String, Option<i32>>>()
}

// Union of bed regions with tolerance.
struct SuperRegion {
    region: Region,
    bed_regions: Vec<Region>,
}

// Creates super regions from multiple bed regions, if they are close enough
// (within tolerance). Later the array of size of this super region will be used
// to calculate depth for each one of the nested bed regions.
fn merge_regions(
    regions: &Vec<(String, u32, u32)>,
    tolerance: u32,
) -> HashMap<String, Vec<SuperRegion>> {
    let ref_id_groups = regions
        .iter()
        .map(|a| (a.0.clone(), (a.1..=a.2)))
        .into_iter()
        .into_group_map();
    let mut ret = HashMap::<String, Vec<SuperRegion>>::new();
    for (ref_id, mut bed_regions) in ref_id_groups.into_iter() {
        bed_regions.sort_by(|a, b| a.start().cmp(b.start()));
        let mut consumed_regions = vec![bed_regions.first().unwrap().clone()];
        let mut super_start = *bed_regions.first().unwrap().start();
        let mut super_end = *bed_regions.first().unwrap().end();
        for range in bed_regions.into_iter().skip(1) {
            if *range.start() > super_end + tolerance {
                ret.entry(ref_id.clone()).or_default().push(SuperRegion {
                    region: super_start..=super_end,
                    bed_regions: consumed_regions,
                });
                consumed_regions = Vec::<Region>::new();
                super_start = *range.start();
            }
            super_end = std::cmp::max(super_end, *range.end());
            consumed_regions.push(range);
        }
        ret.entry(ref_id.clone()).or_default().push(SuperRegion {
            region: super_start..=super_end,
            bed_regions: consumed_regions,
        });
    }
    ret
}

fn get_refid_bounds(
    mut ref_ids: Vec<i32>,
    reader: &mut Reader,
) -> HashMap<i32, Option<Range<usize>>> {
    ref_ids.sort();
    let mut refid_bounds = HashMap::<i32, Option<Range<usize>>>::new();
    let meta_blocks = reader.file_meta.view_blocks(&Fields::RefID);
    // If stats were collected for this file it's possible to narrow binary search.
    if meta_blocks.first().unwrap().stats.is_some() {
        dbg!(meta_blocks.len());
        let mut meta_iter = meta_blocks.iter();
        let mut cur_offset: usize = 0;
        for ref_id in ref_ids.into_iter() {
            // Iterate thorugh blocks, until one with max value bigger than cur_chr is found.
            while let Some(meta_block) = meta_iter.next() {
                if &meta_block.stats.as_ref().unwrap().max_value >= &ref_id {
                    if &meta_block.stats.as_ref().unwrap().min_value > &ref_id {
                        refid_bounds.insert(ref_id, None);
                    } else {
                        // This block may contain the REFID. Later on we will search this block to determine this ultimately.
                        refid_bounds.insert(
                            ref_id,
                            Some(cur_offset..(cur_offset + meta_block.numitems as usize)),
                        );
                    }
                    break;
                }
                cur_offset += meta_block.numitems as usize;
            }
        }
    }
    refid_bounds
}

pub trait DepthWrite {
    fn write_depth(&self, depth: &(i32, u32, u32));
}

fn process_depth_query<W: DepthWrite>(
    reader: &mut Reader,
    output: &mut W,
    buf: &mut Vec<i32>,
    super_regions: Vec<SuperRegion>,
    first_record: usize,
    ref_id: i32,
) {
    let mut counter = 0;
    let mut gbam_buffer_record = GbamRecord::default();
    gbam_buffer_record.cigar = Some(super::cigar::Cigar(Vec::with_capacity(100000)));
    for SuperRegion {
        region: super_region,
        bed_regions,
    } in super_regions.into_iter()
    {
        counter += 1;
        // dbg!(counter);
        buf.clear();
        buf.resize((super_region.end() - super_region.start() + 1) as usize, 0);
        calc_depth(
            reader,
            first_record,
            ref_id,
            &super_region,
            buf,
            &mut gbam_buffer_record,
        );
        let offset = *super_region.start();
        for bed_region in bed_regions.into_iter() {
            for pos in bed_region {
                debug_assert!(pos >= offset);
                output.write_depth(&(ref_id, pos, buf[(pos - offset) as usize] as u32));
            }
        }
    }
}

struct ConsolePrinter {
    ref_id_to_chr: HashMap<i32, String>,
}
impl ConsolePrinter {
    pub fn new(ref_id_to_chr: HashMap<i32, String>) -> Self {
        Self { ref_id_to_chr }
    }
}
impl DepthWrite for ConsolePrinter {
    fn write_depth(&self, depth: &(i32, u32, u32)) {
        println!(
            "{:?}\t{}\t{}",
            self.ref_id_to_chr.get(&depth.0),
            depth.1,
            depth.2
        );
    }
}

pub fn get_regions_depths(reader: &mut Reader, regions: &Vec<(String, u32, u32)>) {
    let ref_id_to_chr = reader
    .file_meta
    .get_ref_seqs()
    .iter().enumerate()
    .map(|(k, (refid, _))| (k as i32, refid.clone()))
    .collect();
    let mut printer = ConsolePrinter::new(ref_id_to_chr);
    get_regions_depths_with_printer(reader, regions, &mut printer);
}

// Approach as in https://github.com/brentp/mosdepth
/// Get depth at position. The file should be sorted!
pub fn get_regions_depths_with_printer<W: DepthWrite>(
    reader: &mut Reader,
    regions: &Vec<(String, u32, u32)>,
    output: &mut W,
) {
    if !reader
        .parsing_template
        .check_if_active(&[Fields::RefID, Fields::Pos, Fields::RawCigar])
    {
        panic!("The reader should have parsing template which includes REFID, POS and CIGAR.");
    }

    let super_regions = merge_regions(regions, 300_000);
    let ref_id_to_chr = get_chr_name_mapping(super_regions.keys(), reader);
    let depth_queries = super_regions
        .into_iter()
        .map(|(k, v)| (ref_id_to_chr.get(&k).unwrap().unwrap(), v))
        .collect::<BTreeMap<i32, Vec<SuperRegion>>>();

    let ref_ids_bounds =
        get_refid_bounds(ref_id_to_chr.values().copied().flatten().collect(), reader);

    let mut buf = Vec::<i32>::new();
    for (ref_id, super_regions) in depth_queries.into_iter() {
        if let Some(ref_id_pos_hint) = ref_ids_bounds.get(&ref_id).unwrap() {
            // Find number of record when records with this ref_id begin.
            // We don't need to fetch anything except refid to find first record with requested refid, so disable other fields fetching.
            reader.fetch_only(&[Fields::RefID]);
            let first_record = find_refid(reader, ref_id, ref_id_pos_hint);
            reader.restore_template();
            if let Some(first_pos) = first_record {
                process_depth_query(reader, output, &mut buf, super_regions, first_pos, ref_id);
            }
        }
    }
}

fn find_leftmost_block(id: i32, block_metas: &Vec<BlockMeta>) -> usize {
    let mut left = 0;
    let mut right = block_metas.len();
    while left < right {
        let mid = (left + right) / 2;
        let max_val = &block_metas[mid].stats.as_ref().unwrap().max_value;
        match max_val.cmp(&id) {
            Ordering::Equal | Ordering::Greater => right = mid,
            Ordering::Less => left = mid + 1,
        }
    }
    left
}

// For each guess there may be I/O operation with decompression, so this method is not fast.
fn find_refid(reader: &mut Reader, chr_num: i32, range: &Range<usize>) -> Option<usize> {
    let pred = |num: usize, buf: &mut GbamRecord| {
        reader.fill_record(num, buf);
        buf.refid.unwrap().cmp(&chr_num)
    };
    // dbg!("Looking for", chr_num);
    binary_search(range.start, range.end, pred)
}

/// Searches for the first record which satisfies predicate.
fn binary_search<F: FnMut(usize, &mut GbamRecord) -> Ordering>(
    mut left: usize,
    mut right: usize,
    mut cmp: F,
) -> Option<usize> {
    let mut buf = GbamRecord::default();
    let end = right;

    while left < right {
        let mid = (left + right) / 2;
        // println!("Mid - {:?}", mid);
        // dbg!(mid, &buf.refid);
        match cmp(mid, &mut buf) {
            Ordering::Less => left = mid + 1,
            Ordering::Equal | Ordering::Greater => right = mid,
        }
    }
    // dbg!("Finished search");
    if left == end || cmp(left, &mut buf) != Ordering::Equal {
        return None;
    }
    Some(left)
}

pub fn calc_depth(
    reader: &mut Reader,
    mut record_num: usize,
    refid: i32,
    super_region: &RangeInclusive<u32>,
    buffer: &mut Vec<i32>,
    buf: &mut GbamRecord,
) {
    let sweeping_line = buffer;
    let mut cur_depth = 0;
    // dbg!("What?");
    loop {
        if record_num >= reader.amount {
            break;
        }

        reader.fill_record(record_num, buf);
        record_num += 1;

        let read_start = buf.pos.unwrap().try_into().unwrap();
        let base_cov = buf.cigar.as_ref().unwrap().base_coverage();
        let read_end = read_start + base_cov;

        // println!("That is test! {} end {}", read_start, read_end);

        if read_start >= *super_region.end() || buf.refid.unwrap() > refid {
            break;
        }

        if super_region.contains(&read_start) {
            sweeping_line[(read_start - super_region.start()) as usize] += 1;
        }
        // Adding 1 since read_end is also covered.
        if super_region.contains(&(read_end + 1)) {
            // Record starts before the region and end on it.
            if !super_region.contains(&read_start) {
                cur_depth += 1;
            }
            sweeping_line[(read_end + 1 - super_region.start()) as usize] -= 1;
        }
    }
    // dbg!("What?");
    let depth_of_region = sweeping_line;

    for sweep in depth_of_region.iter_mut() {
        debug_assert!(!(*sweep < 0 && -*sweep > cur_depth));
        cur_depth += *sweep;
        *sweep = cur_depth;
    }
}

// println!(
//     "REF_ID: {}\nPOS: {}\nbase_cov: {}\ncigar_len: {}\nsearch pos: {}\nEND: {}\nDEPTH: {}\n",
//     buf.refid.unwrap(),
//     read_pos,
//     base_cov,
//     buf.cigar.as_ref().unwrap().0.len(),
//     pos,
//     read_pos + i32::try_from(base_cov).unwrap(),
//     read_pos > pos
// );
