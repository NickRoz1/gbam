use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use byteorder::{LittleEndian, WriteBytesExt};

/// A SAM record consists of fields separated by tabs
#[derive(Debug, Clone)]
pub struct SamRecord {
    pub qname: String,
    pub flag: u16,
    pub rname: String,
    pub pos: u32,
    pub mapq: u8,
    pub cigar: String,
    pub rnext: String,
    pub pnext: u32,
    pub tlen: i32,
    pub seq: String,
    pub qual: String,
    pub tags: Vec<String>,
}

impl SamRecord {
    pub fn from_line(line: &str) -> Option<Self> {
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 11 {
            return None;
        }

        Some(SamRecord {
            qname: parts[0].to_string(),
            flag: parts[1].parse().ok()?,
            rname: parts[2].to_string(),
            pos: parts[3].parse().ok()?,
            mapq: parts[4].parse().ok()?,
            cigar: parts[5].to_string(),
            rnext: parts[6].to_string(),
            pnext: parts[7].parse().ok()?,
            tlen: parts[8].parse().ok()?,
            seq: parts[9].to_string(),
            qual: parts[10].to_string(),
            tags: parts[11..].iter().map(|s| s.to_string()).collect(),
        })
    }
}

#[derive(Debug, Clone)]
pub struct SamHeader {
    pub lines: Vec<String>,
}

impl SamHeader {
    pub fn as_bytes_and_offset(&self) -> io::Result<(Vec<u8>, usize)> {
        let joined = self.lines.join("\n") + "\n";
        let text_bytes = joined.as_bytes();
        let l_text = text_bytes.len();

        let mut bytes = Vec::new();
        {
            // Use a temporary scope and reference to avoid move issues
            let bytes_ref = &mut bytes;
            bytes_ref.write_u32::<LittleEndian>(l_text as u32)?;
            bytes_ref.write_all(text_bytes)?;
        }

        let offset = bytes.len();
        Ok((bytes, offset))
    }

    pub fn reference_names(&self) -> Vec<String> {
        self.lines
            .iter()
            .filter(|line| line.starts_with("@SQ"))
            .filter_map(|line| {
                for field in line.split('\t') {
                    if let Some(name) = field.strip_prefix("SN:") {
                        return Some(name.to_string());
                    }
                }
                None
            })
            .collect()
    }

    pub fn parse_reference_sequences(&self) -> Vec<(String, u32)> {
        self.lines
            .iter()
            .filter(|line| line.starts_with("@SQ"))
            .filter_map(|line| {
                let mut name: Option<String> = None;
                let mut length: Option<u32> = None;

                for field in line.split('\t').skip(1) {
                    if let Some(stripped) = field.strip_prefix("SN:") {
                        name = Some(stripped.to_string());
                    } else if let Some(stripped) = field.strip_prefix("LN:") {
                        length = stripped.parse::<u32>().ok();
                    }
                }

                match (name, length) {
                    (Some(n), Some(l)) => Some((n, l)),
                    _ => None,
                }
            })
            .collect()
    }
}

pub struct SamReader<R: BufRead> {
    reader: R,
    header: SamHeader,
}

impl<R: BufRead> SamReader<R> {
    pub fn new(mut reader: R) -> Self {
        let mut header_lines = Vec::new();
        let mut buffer = String::new();
    
        // Read and collect header lines without consuming the first record line
        let mut position = reader.fill_buf().unwrap();
        while !position.is_empty() {
            let line_end = position.iter().position(|&b| b == b'\n').map(|i| i + 1);
            if let Some(end) = line_end {
                let line = std::str::from_utf8(&position[..end]).unwrap();
                if line.starts_with('@') {
                    header_lines.push(line.trim_end().to_string());
                    reader.consume(end);
                    position = reader.fill_buf().unwrap();
                } else {
                    break;
                }
            } else {
                break;
            }
        }
    
        SamReader {
            reader,
            header: SamHeader { lines: header_lines },
        }
    }

    pub fn header(&self) -> io::Result<(Vec<u8>, usize)> {
        self.header.as_bytes_and_offset()
    }

    pub fn reference_sequences(&self) -> Vec<(String, u32)> {
        self.header.parse_reference_sequences()
    }

    pub fn reference_names(&self) -> Vec<String> {
        self.header.reference_names()
    }

    pub fn records(&mut self) -> SamParsedRecords<'_, R> {
        SamParsedRecords::new(&mut self.reader)
    }

    pub fn raw_records(&mut self) -> SamRawRecords<'_, R> {
        SamRawRecords::new(&mut self.reader)
    }
}

// Yields parsed SamRecord from reader
pub struct SamParsedRecords<'a, R: BufRead> {
    reader: &'a mut R,
    buffer: String,
}

impl<'a, R: BufRead> SamParsedRecords<'a, R> {
    pub fn new(reader: &'a mut R) -> Self {
        Self {
            reader,
            buffer: String::new(),
        }
    }
}

impl<'a, R: BufRead> Iterator for SamParsedRecords<'a, R> {
    type Item = io::Result<SamRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        self.buffer.clear();
        match self.reader.read_line(&mut self.buffer) {
            Ok(0) => None,
            Ok(_) => {
                if self.buffer.starts_with('@') {
                    return self.next(); // Skip header lines
                }
                match SamRecord::from_line(self.buffer.trim_end()) {
                    Some(rec) => Some(Ok(rec)),
                    None => Some(Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "Invalid SAM line",
                    ))),
                }
            }
            Err(e) => Some(Err(e)),
        }
    }
}

// Raw SAM record reader for &[u8] lines
pub struct SamRawRecords<'a, R: BufRead> {
    reader: &'a mut R,
    buffer: Vec<u8>,
    scratch: String,
}

impl<'a, R: BufRead> SamRawRecords<'a, R> {
    pub fn new(reader: &'a mut R) -> Self {
        Self {
            reader,
            buffer: Vec::new(),
            scratch: String::new(),
        }
    }

    pub fn next_rec(&mut self) -> Option<io::Result<&[u8]>> {
        self.buffer.clear();
        self.scratch.clear();

        loop {
            match self.reader.read_line(&mut self.scratch) {
                Ok(0) => return None,
                Ok(_) => {
                    if self.scratch.starts_with('@') {
                        self.scratch.clear();
                        continue;
                    }
                    self.buffer.extend_from_slice(self.scratch.as_bytes());
                    return Some(Ok(&self.buffer));
                }
                Err(e) => return Some(Err(e)),
            }
        }
    }
}

// Convenience function to open a SAM file and return header and parsed record iterator
// pub fn read_sam_file(
//     path: &str,
// ) -> io::Result<((Vec<u8>, usize), SamParsedRecords<BufReader<File>>)> {
//     let file = File::open(path)?;
//     let mut reader = BufReader::new(file);

//     let mut header_lines = Vec::new();
//     let mut buffer = String::new();

//     while let Ok(n) = reader.read_line(&mut buffer) {
//         if n == 0 {
//             break;
//         }
//         if buffer.starts_with('@') {
//             header_lines.push(buffer.trim_end().to_string());
//             buffer.clear();
//         } else {
//             break;
//         }
//     }

//     let header = SamHeader { lines: header_lines };
//     let (raw_header, offset) = header.as_bytes_and_offset();

//     Ok(((raw_header, offset), SamParsedRecords::new(&mut reader)))
// }
