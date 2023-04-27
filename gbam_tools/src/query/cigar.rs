use std::{fmt::Display, slice::Iter};

use bam_tools::record::bamrawrecord::decode_cigar;
use byteorder::ByteOrder;
use byteorder::LittleEndian;
use byteorder::WriteBytesExt;
use std::io::Write;

#[derive(Debug, Clone)]
pub struct Op(pub u32);

impl Op {
    pub fn new(val: u32) -> Op {
        Op(val)
    }
    /// True if operation is one of M, =, X, D, N
    pub fn is_consuming_reference(&self) -> bool {
        let op = self.0 & 0xF;
        match op {
            0 | 2 | 3 | 7 | 8 => true,
            _ => false,
        }
    }
    /// Length of operator
    pub fn length(&self) -> u32 {
        return self.0 >> 4;
    }

    /// Type of operator itself
    pub fn op_type(&self) -> char {
        match self.0 & 0xF {
            0 => 'M',
            1 => 'I',
            2 => 'D',
            3 => 'N',
            4 => 'S',
            5 => 'H',
            6 => 'P',
            7 => '=',
            8 => 'X',
            _ => panic!("Unexpected cigar operation"),
        }
    }
}

#[derive(Debug)]
pub struct Cigar(pub Vec<Op>);

impl Cigar {
    pub fn new(ops: Vec<Op>) -> Cigar {
        Cigar(ops)
    }

    pub fn base_coverage(&self) -> u32 {
        let mut count = 0;
        for op in self.ops() {
            if op.is_consuming_reference() {
                count += op.length();
            }
        }
        count
    }

    pub fn ops(&self) -> Iter<Op> {
        self.0.iter()
    }

    pub fn write_as_bytes<T: ByteOrder>(&self, bytes: &mut Vec<u8>) {
        self.ops()
            .for_each(|op| bytes.write_u32::<T>(op.0).unwrap());
    }
}

// impl Display for Cigar {
//     fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
//         write!(f, "{}", decode_cigar(&self.0[..]))
//     }
// }
