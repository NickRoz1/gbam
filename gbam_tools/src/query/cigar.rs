use std::{slice::Iter};

use serde::{Serialize, Deserialize};


use byteorder::ByteOrder;
use byteorder::WriteBytesExt;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Op(pub u32);

impl Op {
    pub fn new(val: u32) -> Op {
        Op(val)
    }
    /// True if operation is one of M, =, X, D, N
    pub fn is_consuming_reference(&self) -> bool {
        let op = self.0 & 0xF;
        matches!(op, 0 | 2 | 3 | 7 | 8)
    }
    /// Returns whether the operation kind causes the alignment to consume the read.
    pub fn consumes_read(&self) -> bool {
        let op = self.0 & 0xF;
        matches!(op, 0 | 1 | 4 | 7 | 8)
    }
    /// Length of operator
    pub fn length(&self) -> u32 {
        self.0 >> 4
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

#[derive(Debug, Serialize, Deserialize)]
pub struct Cigar(pub Vec<Op>);

pub fn base_coverage(arr: &[Op]) -> u32 {
    let mut count = 0;
    for op in arr {
        if op.is_consuming_reference() {
            count += op.length();
        }
    }
    count
}

impl Cigar {
    pub fn new(ops: Vec<Op>) -> Cigar {
        Cigar(ops)
    }

    /// Calculates the read length.
    ///
    /// This sums the lengths of the CIGAR operations that consume the read, i.e., alignment
    /// matches (`M`), insertions to the reference (`I`), soft clips (`S`), sequence matches (`=`),
    /// and sequence mismatches (`X`).
    pub fn read_length(&self) -> u32 {
        let mut count = 0;
        for op in self.ops() {
            if op.consumes_read() {
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

impl std::fmt::Display for Cigar {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.ops().try_for_each(|op| {
            write!(f, "{}{}", op.length(), op.op_type())
        })
    }
}
