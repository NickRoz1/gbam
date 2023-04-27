use crate::reader::parse_tmplt::ParsingTemplate;
use crate::reader::records::Records;
use rust_htslib::bam;
use std::io::Write;

use std::convert::TryFrom;
use std::fs::File;

/// Converts GBAM file to BAM file. This uses the `noodles bam writer`.
pub fn gbam_to_bam(in_path: &str, out_path: &str) {
    let file = File::open(in_path).unwrap();
    let mut template = ParsingTemplate::new();
    template.set_all();
    let mut reader = crate::reader::reader::Reader::new(file, template).unwrap();

    let mut bam_header = bam::Header::new();
    let ref_seqs = reader.file_meta.get_ref_seqs();

    for ref_seq in ref_seqs {
        bam_header.push_record(
            bam::header::HeaderRecord::new(b"SQ")
                .push_tag(b"SN", &ref_seq.0)
                .push_tag(b"LN", &ref_seq.1),
        );
    }

    let mut records_it = Records::new(&mut reader);

    let mut out = bam::Writer::from_path(out_path, &bam_header, bam::Format::Bam).unwrap();
    out.set_threads(4).unwrap();

    let mut cigar_buf = Vec::new();
    while let Some(rec) = records_it.next_rec() {
        let mut record = bam::Record::new();

        record.set_bin(rec.bin.unwrap());
        record.set_tid(rec.refid.unwrap());
        record.set_mapq(rec.mapq.unwrap());
        record.set_pos(rec.pos.unwrap() as i64);
        record.set_flags(rec.flag.unwrap());
        record.set_mtid(rec.next_ref_id.unwrap());
        record.set_mpos(rec.next_pos.unwrap() as i64);
        record.set_insert_size(rec.tlen.unwrap() as i64);
        let rec_seq_len = rec.seq.as_ref().unwrap().len();
        let mut qual = rec.qual.as_ref().unwrap().clone();
        if qual.is_empty() {
            qual = vec![255; rec_seq_len];
        }

        cigar_buf.clear();
        rec.cigar.as_ref().unwrap().ops().for_each(|op| {
            cigar_buf
                .write_all(op.length().to_string().as_bytes())
                .unwrap();
            cigar_buf.push(op.op_type() as u8);
        });

        let bam_cigar = bam::record::CigarString::try_from(&cigar_buf[..]).unwrap();
        record.set_data(&rec.tags.as_ref().unwrap()[..]);
        record.set(
            &rec.read_name.as_ref().unwrap()[..rec.read_name.as_ref().unwrap().len() - 1],
            Some(&bam_cigar),
            rec.seq.as_ref().unwrap().as_bytes(),
            &qual[..],
        );

        out.write(&record).unwrap();
    }
}
