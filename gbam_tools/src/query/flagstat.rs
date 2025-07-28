use crate::reader::parse_tmplt::ParsingTemplate;
use crate::reader::reader::Reader;
use crate::reader::record::GbamRecord;
use bam_tools::record::fields::Fields;
use bitflags::bitflags;
use rayon::prelude::*;
use std::fmt;
use std::fs::File;
use std::str;
use std::string::String;

// https://github.com/samtools/htslib/blob/32de287eafdafc45dde0a22244b72697294f161d/htslib/sam.h
bitflags! {
    #[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
        struct BamFlags: u32 {
            /// @abstract the read is paired in sequencing, no matter whether it is mapped in a pair.
            const BAM_FPAIRED =        1;
            /// @abstract the read is mapped in a proper pair.
            const BAM_FPROPER_PAIR =   2;
            /// @abstract the read itself is unmapped; conflictive with BAM_FPROPER_PAIR.
            const BAM_FUNMAP =         4;
            /// @abstract the mate is unmapped.
            const BAM_FMUNMAP =        8;
            /// @abstract the read is mapped to the reverse strand.
            const BAM_FREVERSE =      16;
            /// @abstract the mate is mapped to the reverse strand.
            const BAM_FMREVERSE =     32;
            /// @abstract this is read1.
            const BAM_FREAD1 =        64;
            /// @abstract this is read2.
            const BAM_FREAD2 =       128;
            /// @abstract not primary alignment.
            const BAM_FSECONDARY =   256;
            /// @abstract QC failure.
            const BAM_FQCFAIL =      512;
            /// @abstract optical or PCR duplicate.
            const BAM_FDUP =        1024;
            /// @abstract supplementary alignment.
            const BAM_FSUPPLEMENTARY = 2048;
    }
}

#[derive(Default)]
struct Stats {
    pub n_reads: [i64; 2],
    pub n_mapped: [i64; 2],
    pub n_pair_all: [i64; 2],
    pub n_pair_map: [i64; 2],
    pub n_pair_good: [i64; 2],
    pub n_sgltn: [i64; 2],
    pub n_read1: [i64; 2],
    pub n_read2: [i64; 2],
    pub n_dup: [i64; 2],
    pub n_diffchr: [i64; 2],
    pub n_diffhigh: [i64; 2],
    pub n_secondary: [i64; 2],
    pub n_supp: [i64; 2],
    pub n_primary: [i64; 2],
    pub n_pmapped: [i64; 2],
    pub n_pdup: [i64; 2],
}

impl Stats {
    fn add_two_arrs(dest: &mut [i64; 2], src: &[i64; 2]) {
        dest[0] += src[0];
        dest[1] += src[1];
    }
    pub fn add(&mut self, other: &Stats) {
        Self::add_two_arrs(&mut self.n_reads, &other.n_reads);
        Self::add_two_arrs(&mut self.n_mapped, &other.n_mapped);
        Self::add_two_arrs(&mut self.n_pair_all, &other.n_pair_all);
        Self::add_two_arrs(&mut self.n_pair_map, &other.n_pair_map);
        Self::add_two_arrs(&mut self.n_pair_good, &other.n_pair_good);
        Self::add_two_arrs(&mut self.n_sgltn, &other.n_sgltn);
        Self::add_two_arrs(&mut self.n_read1, &other.n_read1);
        Self::add_two_arrs(&mut self.n_read2, &other.n_read2);
        Self::add_two_arrs(&mut self.n_dup, &other.n_dup);
        Self::add_two_arrs(&mut self.n_diffchr, &other.n_diffchr);
        Self::add_two_arrs(&mut self.n_diffhigh, &other.n_diffhigh);
        Self::add_two_arrs(&mut self.n_secondary, &other.n_secondary);
        Self::add_two_arrs(&mut self.n_supp, &other.n_supp);
        Self::add_two_arrs(&mut self.n_primary, &other.n_primary);
        Self::add_two_arrs(&mut self.n_pmapped, &other.n_pmapped);
        Self::add_two_arrs(&mut self.n_pdup, &other.n_pdup);
    }
}

fn percent(n: i64, total: i64) -> String {
    if total != 0 {
        format!("{:.2}%", (n as f64) / (total as f64) * 100.0)
    } else {
        String::from("N/A")
    }
}

impl fmt::Display for Stats {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(
            f,
            "{} + {} in total (QC-passed reads + QC-failed reads)",
            self.n_reads[0], self.n_reads[1]
        )
        .unwrap();
        writeln!(f, "{} + {} primary", self.n_primary[0], self.n_primary[1]).unwrap();
        writeln!(
            f,
            "{} + {} secondary",
            self.n_secondary[0], self.n_secondary[1]
        )
        .unwrap();
        writeln!(f, "{} + {} supplementary", self.n_supp[0], self.n_supp[1]).unwrap();
        writeln!(f, "{} + {} duplicates", self.n_dup[0], self.n_dup[1]).unwrap();
        writeln!(
            f,
            "{} + {} primary duplicates",
            self.n_pdup[0], self.n_pdup[1]
        )
        .unwrap();
        writeln!(
            f,
            "{} + {} mapped ({} : {})",
            self.n_mapped[0],
            self.n_mapped[1],
            percent(self.n_mapped[0], self.n_reads[0]),
            percent(self.n_mapped[1], self.n_reads[1])
        )
        .unwrap();
        writeln!(
            f,
            "{} + {} primary mapped ({} : {})",
            self.n_pmapped[0],
            self.n_pmapped[1],
            percent(self.n_pmapped[0], self.n_primary[0]),
            percent(self.n_pmapped[1], self.n_primary[1])
        )
        .unwrap();
        writeln!(
            f,
            "{} + {} paired in sequencing",
            self.n_pair_all[0], self.n_pair_all[1]
        )
        .unwrap();
        writeln!(f, "{} + {} read1", self.n_read1[0], self.n_read1[1]).unwrap();
        writeln!(f, "{} + {} read2", self.n_read2[0], self.n_read2[1]).unwrap();
        writeln!(
            f,
            "{} + {} properly paired ({} : {})",
            self.n_pair_good[0],
            self.n_pair_good[1],
            percent(self.n_pair_good[0], self.n_pair_all[0]),
            percent(self.n_pair_good[1], self.n_pair_all[1])
        )
        .unwrap();
        writeln!(
            f,
            "{} + {} with itself and mate mapped",
            self.n_pair_map[0], self.n_pair_map[1]
        )
        .unwrap();
        writeln!(
            f,
            "{} + {} singletons ({} : {})",
            self.n_sgltn[0],
            self.n_sgltn[1],
            percent(self.n_sgltn[0], self.n_pair_all[0]),
            percent(self.n_sgltn[1], self.n_pair_all[1])
        )
        .unwrap();
        writeln!(
            f,
            "{} + {} with mate mapped to a different chr",
            self.n_diffchr[0], self.n_diffchr[1]
        )
        .unwrap();
        write!(
            f,
            "{} + {} with mate mapped to a different chr (mapQ>=5)",
            self.n_diffhigh[0], self.n_diffhigh[1]
        )
    }
}

fn collect(rec: &GbamRecord, stats: &mut Stats) {
    let record_flag = BamFlags::from_bits(rec.flag.unwrap() as u32).unwrap();
    let w = record_flag.contains(BamFlags::BAM_FQCFAIL) as usize;

    stats.n_reads[w] += 1;

    if record_flag.contains(BamFlags::BAM_FSECONDARY) {
        stats.n_secondary[w] += 1;
    } else if record_flag.contains(BamFlags::BAM_FSUPPLEMENTARY) {
        stats.n_supp[w] += 1;
    } else {
        stats.n_primary[w] += 1;

        if record_flag.contains(BamFlags::BAM_FPAIRED) {
            stats.n_pair_all[w] += 1;
            if record_flag.contains(BamFlags::BAM_FPROPER_PAIR)
                && !record_flag.contains(BamFlags::BAM_FUNMAP)
            {
                stats.n_pair_good[w] += 1;
            }
            if record_flag.contains(BamFlags::BAM_FREAD1) {
                stats.n_read1[w] += 1;
            }
            if record_flag.contains(BamFlags::BAM_FREAD2) {
                stats.n_read2[w] += 1;
            }
            if record_flag.contains(BamFlags::BAM_FMUNMAP)
                && !record_flag.contains(BamFlags::BAM_FUNMAP)
            {
                stats.n_sgltn[w] += 1;
            }
            if !record_flag.contains(BamFlags::BAM_FUNMAP)
                && !record_flag.contains(BamFlags::BAM_FMUNMAP)
            {
                stats.n_pair_map[w] += 1;
                if rec.next_ref_id.unwrap() != rec.refid.unwrap() {
                    stats.n_diffchr[w] += 1;
                    if rec.mapq.unwrap() >= 5 {
                        stats.n_diffhigh[w] += 1;
                    }
                }
            }
        }
        if !record_flag.contains(BamFlags::BAM_FUNMAP) {
            stats.n_pmapped[w] += 1;
        }
        if record_flag.contains(BamFlags::BAM_FDUP) {
            stats.n_pdup[w] += 1;
        }
    }
    if !record_flag.contains(BamFlags::BAM_FUNMAP) {
        stats.n_mapped[w] += 1;
    }
    if record_flag.contains(BamFlags::BAM_FDUP) {
        stats.n_dup[w] += 1;
    }
}

pub fn collect_stats(file: File) {
    let tmplt = ParsingTemplate::new();
    let reader = Reader::new(file.try_clone().unwrap(), tmplt).unwrap();
    let total_records = reader.amount;
    let file_meta = reader.file_meta;

    let file_stats = (0..total_records)
        .into_par_iter()
        .chunks(500_000)
        .map(|records_range| {
            let mut stats = Stats::default();

            let mut rec = GbamRecord::default();
            let mut tmplt = ParsingTemplate::new();
            tmplt.set(&Fields::Flags, true);
            tmplt.set(&Fields::RefID, true);
            tmplt.set(&Fields::NextRefID, true);
            tmplt.set(&Fields::Mapq, true);

            let mut reader =
                Reader::new_with_meta(file.try_clone().unwrap(), tmplt, &file_meta, None).unwrap();

            for rec_num in records_range {
                reader.fill_record(rec_num, &mut rec);
                collect(&rec, &mut stats);
            }

            stats
        })
        .reduce(Stats::default, |mut a, b| {
            a.add(&b);
            a
        });

    println!("{file_stats}");
}
