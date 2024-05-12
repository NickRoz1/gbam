#[allow(dead_code)]
// RFC 1952 § 2.3.1
pub(crate) const MAGIC_NUMBER: [u8; 2] = [0x1f, 0x8b];

#[allow(dead_code)]
pub(crate) const MTIME_NONE: u32 = 0;

// ID1 (1) + ID2 (1) + CM (1) + FLG (1) + MTIME (4) + XLF (1) + OS (1)
pub(crate) const HEADER_SIZE: usize = 10;

// CRC32 (4) + ISIZE (4)
pub(crate) const TRAILER_SIZE: usize = 8;

// XLEN (2)
const GZIP_XLEN_SIZE: usize = 2;

// SI1 (1) + SI2 (1) + SLEN (2) + BSIZE (2)
const BGZF_XLEN: usize = 6;

pub(crate) const BGZF_HEADER_SIZE: usize = HEADER_SIZE + GZIP_XLEN_SIZE + BGZF_XLEN;

#[non_exhaustive]
#[allow(dead_code)]
pub(crate) enum CompressionMethod {
    Deflate = 8,
}

#[non_exhaustive]
#[allow(dead_code)]
pub(crate) enum OperatingSystem {
    Unknown = 255,
}
