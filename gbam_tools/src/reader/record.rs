use std::io::Write;

use itertools::Itertools;
use serde::{Serialize, Deserialize};

use bam_tools::record::{
    bamrawrecord::{decode_seq, put_sequence},
    fields::Fields,
};

use crate::query::cigar::base_coverage;
use byteorder::{LittleEndian, ReadBytesExt, WriteBytesExt};
use std::mem;


use crate::{query::cigar::Cigar, query::cigar::Op, U32_SIZE};


#[derive(Debug, Default, Serialize, Deserialize)]
/// Represents a GBAM record in which some fields may be omitted.
pub struct GbamRecord {
    /// Reference sequence ID
    pub refid: Option<i32>,
    /// 0-based leftmost coordinate
    pub pos: Option<i32>,
    /// Mapping quality
    pub mapq: Option<u8>,
    /// BAI index bin,
    pub bin: Option<u16>,
    /// Bitwise flags
    pub flag: Option<u16>,
    /// Ref-ID of the next segment
    pub next_ref_id: Option<i32>,
    /// 0-based leftmost pos of the next segmen
    pub next_pos: Option<i32>,
    /// Template length
    pub tlen: Option<i32>,
    /// Read name
    pub read_name: Option<Vec<u8>>,
    /// CIGAR
    pub cigar: Option<Cigar>,
    /// 4-bit  encoded  read
    pub seq: Option<String>,
    /// Phred-scaled base qualities.
    pub qual: Option<Vec<u8>>,
    /// List of auxiliary data
    pub tags: Option<Vec<u8>>,
}

pub fn parse_cigar(bytes: &[u8], prealloc: &mut Cigar) {
    prealloc.0.resize(bytes.len() / U32_SIZE, Op::new(0));
    for (i, mut chunk) in bytes.chunks(U32_SIZE).enumerate() {
        prealloc.0[i] = Op::new(chunk.read_u32::<LittleEndian>().unwrap());
    }
}

// TODO :: ADD TEMPLATE LENGTHS TO GBAM RECORD
// TODO :: REMOVE CG TAG FROM ORIGINAL FILE
impl GbamRecord {
    pub(crate) fn parse_from_bytes(&mut self, field: &Fields, mut bytes: &[u8]) {
        match field {
            Fields::RefID => self.refid = Some(bytes.read_i32::<LittleEndian>().unwrap()),
            Fields::Pos => self.pos = Some(bytes.read_i32::<LittleEndian>().unwrap()),
            Fields::Mapq => self.mapq = Some(bytes[0].to_owned()),
            Fields::Bin => self.bin = Some(bytes.read_u16::<LittleEndian>().unwrap()),
            Fields::Flags => self.flag = Some(bytes.read_u16::<LittleEndian>().unwrap()),
            Fields::NextRefID => self.next_ref_id = Some(bytes.read_i32::<LittleEndian>().unwrap()),
            Fields::NextPos => self.next_pos = Some(bytes.read_i32::<LittleEndian>().unwrap()),
            Fields::TemplateLength => self.tlen = Some(bytes.read_i32::<LittleEndian>().unwrap()),
            Fields::ReadName => self.read_name = Some(bytes.to_vec()),
            Fields::RawCigar => {
                parse_cigar(bytes, self.cigar.get_or_insert(Cigar::new(Vec::new())));
            }
            Fields::RawSequence => {
                decode_seq(bytes, self.seq.get_or_insert(String::new()))
            },
            Fields::RawQual => self.qual = Some(bytes.to_vec()),
            Fields::RawTags => self.tags = Some(bytes.to_vec()),
            _ => panic!("Not yet covered type: {}", field),
        }
    }

    /// Only support full records. Do not call if the GBAM record is not fully filled.
    ///
    /// Layout:
    ///
    /// block_size                       uint32_t
    /// refID                            int32_t
    /// pos                              int32 t
    /// l_read_name                      uint8_t
    /// mapq                             uint8_t
    /// bin                              uint16_t
    /// n cigar op                       uint16_t
    /// flag                             uint16_t
    /// l_seq                            uint32_t
    /// next_refid                       int32_t
    /// next pos                         int32_t
    /// tlen                             int32_t
    /// read name                        char[l_read_name]
    /// cigar                            uint32_t[n_cigar_op]
    /// seq                              uint8_t[(l_seq+1)/2]
    /// qual                             char[l_seq]
    /// UNTIL THE END OF BLOCK:
    /// tag                              char[2]
    /// val_type                         char
    /// tag_value                        by_val_type
    pub fn convert_to_bytes(&self, bytes: &mut Vec<u8>) {
        let n_byte = mem::size_of::<u32>()
            + mem::size_of::<i32>() * 2
            + mem::size_of::<u8>() * 2
            + mem::size_of::<u16>() * 3
            + mem::size_of::<u32>()
            + mem::size_of::<i32>() * 3
            + self.cigar.as_ref().unwrap().0.len() * mem::size_of::<u32>()
            + self.read_name.as_ref().unwrap().len()
            + (self.seq.as_ref().unwrap_or(&String::new()).len()+1)/2
            + self.qual.as_ref().unwrap_or(&Vec::new()).len()
            + self.tags.as_ref().unwrap().len();

        bytes.resize(n_byte, 0);

        (&mut bytes[0..4])
            .write_u32::<LittleEndian>((n_byte - mem::size_of::<u32>()) as u32)
            .unwrap();
        (&mut bytes[4..8])
            .write_i32::<LittleEndian>(self.refid.unwrap())
            .unwrap();
        (&mut bytes[8..12])
            .write_i32::<LittleEndian>(self.pos.unwrap())
            .unwrap();
        (&mut bytes[12..13])
            .write_u8(self.read_name.as_ref().unwrap().len() as u8)
            .unwrap();
        (&mut bytes[13..14]).write_u8(self.mapq.unwrap()).unwrap();
        (&mut bytes[14..16])
            .write_u16::<LittleEndian>(self.bin.unwrap())
            .unwrap();
        (&mut bytes[16..18])
            .write_u16::<LittleEndian>(self.cigar.as_ref().unwrap().0.len() as u16)
            .unwrap();
        (&mut bytes[18..20])
            .write_u16::<LittleEndian>(self.flag.unwrap())
            .unwrap();
        // Since we can't get right value from SEQ field (look in SAM/BAM documentation).;
        (&mut bytes[20..24])
            .write_u32::<LittleEndian>(self.qual.as_ref().unwrap_or(&Vec::new()).len() as u32)
            .unwrap();
        (&mut bytes[24..28])
            .write_i32::<LittleEndian>(self.next_ref_id.unwrap())
            .unwrap();
        (&mut bytes[28..32])
            .write_i32::<LittleEndian>(self.next_pos.unwrap())
            .unwrap();
        (&mut bytes[32..36])
            .write_i32::<LittleEndian>(self.tlen.unwrap())
            .unwrap();
        let unsized_data = &mut bytes[36..];
        let (mut read_name, unsized_data) =
            unsized_data.split_at_mut(self.read_name.as_ref().unwrap().len());
        read_name
            .write_all(self.read_name.as_ref().unwrap())
            .unwrap();
        let (cigar, unsized_data) =
            unsized_data.split_at_mut(self.cigar.as_ref().unwrap().0.len() * mem::size_of::<u32>());
        self.cigar
            .as_ref()
            .unwrap()
            .ops()
            .zip_eq(cigar.chunks_mut(mem::size_of::<u32>()))
            .for_each(|(op, mut buf)| buf.write_u32::<LittleEndian>(op.0).unwrap());
        let seq_len = (self.seq.as_ref().unwrap_or(&String::new()).len()+1)/2;
        let (seq, unsized_data) = unsized_data.split_at_mut(seq_len);
        put_sequence(seq, self.seq.as_ref().unwrap_or(&String::new()).len(), self.seq.as_ref().unwrap_or(&String::new())).unwrap();
        let (mut qual, mut unsized_data) =
            unsized_data.split_at_mut(self.qual.as_ref().unwrap_or(&Vec::new()).len());
        qual.write_all(self.qual.as_ref().unwrap_or(&Vec::new())).unwrap();
        assert!(unsized_data.len() == self.tags.as_ref().unwrap().len());
        unsized_data
            .write_all(self.tags.as_ref().unwrap())
            .unwrap();
        assert!(unsized_data.is_empty());
    }

    /// Write tags into a byte buffer.
    pub fn convert_tags_to_bytes(&self, bytes: &mut Vec<u8>) {
        let n_byte = self.tags.as_ref().unwrap().len();

        bytes.reserve(n_byte + bytes.len());
        bytes.write_all(self.tags.as_ref().unwrap()).unwrap();
    }

    /// Returns the alignment span.
    pub fn alignment_span(&self) -> u32 {
        base_coverage(&self.cigar.as_ref().unwrap().0[..])
    }

    /// Returns the alignment start.
    pub fn alignment_start(&self) -> Option<u32> {
        Option::from(self.pos.unwrap() as u32)
    }

    /// Calculates the end position.
    pub fn alignment_end(&self) -> Option<u32> {
        self.alignment_start().and_then(|alignment_start| {
            Option::from(alignment_start + self.alignment_span() - 1)
        })
    }

    pub fn is_reverse(&self) -> bool {
        let flag = self.flag.unwrap();
        (flag & 0x10) == 0x10 as u16
    }

    pub fn is_reverse_complemented(&self) -> bool {
        let flag = self.flag.unwrap();
        (flag & rust_htslib::htslib::BAM_FREVERSE as u16) == rust_htslib::htslib::BAM_FREVERSE as u16
    }

    pub fn is_unmapped(&self) -> bool {
        let flag = self.flag.unwrap();
        (flag & rust_htslib::htslib::BAM_FUNMAP as u16) == rust_htslib::htslib::BAM_FUNMAP as u16
    }
}

impl std::fmt::Display for GbamRecord {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        std::fmt::Debug::fmt(self, f)
    }
}