use std::fs::File;
use std::io::{self, BufRead, Read};
use std::path::Path;
use std::{collections::HashMap};

const SEPARATOR: char = '\t';

/// Source: https://github.com/zaeleus/noodles/blob/90e70874eaa6dd41ac8339933d6dd95bd98080c2/noodles-tabix/examples/tabix_write.rs#L24
fn parse_record(s: &str) -> io::Result<(String, u32, u32)> {
    let mut components = s.splitn(3, SEPARATOR);

    let reference_sequence_name = components
        .next()
        .ok_or_else(|| io::Error::from(io::ErrorKind::InvalidData))?
        .to_owned();

    let start = components
        .next()
        .and_then(|t| t.parse().ok())
        .ok_or_else(|| io::Error::from(io::ErrorKind::InvalidData))?;

    let end = components
        .next()
        .and_then(|t| t.parse().ok())
        .ok_or_else(|| io::Error::from(io::ErrorKind::InvalidData))?;

    if end < start {
        return Err(io::Error::from(io::ErrorKind::InvalidData));
    }

    Ok((reference_sequence_name, start, end))
}

pub fn parse_bed_from_file(path: &Path) -> io::Result<HashMap::<String, Vec<(u32, u32)>>> {
    let mut file = File::open(path)?;
    parse_bed(&mut file)
}

pub fn parse_bed<R: Read>(source: &mut R) -> io::Result<HashMap::<String, Vec<(u32, u32)>>> {
    let mut res = HashMap::<String, Vec<(u32, u32)>>::new();
    let lines = read_lines(source)?;

    for line in lines {
        if let Ok(rec) = parse_record(&line?) {
            res.entry(rec.0).or_insert(Vec::new()).push((rec.1, rec.2));
        }
    }

    Ok(res)
}

fn read_lines<R>(source: &mut R) -> io::Result<io::Lines<io::BufReader<&mut R>>>
where
    R: Read,
{
    Ok(io::BufReader::new(source).lines())
}

pub fn parse_region_query(q: &str) -> io::Result<(&str, u32, u32)> {
    let mut parts = q.split(':');

    let ref_id = parts.next().unwrap();

    let mut region_parts = parts.next().unwrap().split('-');

    let left = region_parts
        .next()
        .unwrap()
        .parse::<u32>()
        .ok()
        .ok_or_else(|| io::Error::from(io::ErrorKind::InvalidData))?;
    let right = region_parts
        .next()
        .unwrap()
        .parse::<u32>()
        .ok()
        .ok_or_else(|| io::Error::from(io::ErrorKind::InvalidData))?;

    if region_parts.next().is_some() || right < left {
        return Err(io::Error::from(io::ErrorKind::InvalidData));
    };

    Ok((ref_id, left, right))
}

pub fn parse_region_query_owned(q: &str) -> io::Result<(String, u32, u32)> {
    let res = parse_region_query(q)?;
    Ok((res.0.to_owned(), res.1, res.2))
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;

    #[test]
    fn test_region_query() {
        let res_1 = parse_region_query("chr3:51289-19238568").unwrap();
        let res_2 = parse_region_query("chrX:346798-23689090").unwrap();
        let res_3 = parse_region_query("chrX:346798-23689090-");
        let res_4 = parse_region_query("chrX:346798-1");
        assert_eq!(res_1, ("chr3", 51289, 19238568));
        assert_eq!(res_2, ("chrX", 346798, 23689090));
        assert!(res_3.is_err());
        assert!(res_4.is_err());
    }

    #[test]
    fn test_bed_reader() {
        let source = "chr3\t51289\t19238568\n\
        chrX\t346798\t23689090-\n\
        chrX\t346798\t23689090
        ";
        let mut reader = Cursor::new(source);
        let res = parse_bed(&mut reader).unwrap();
        assert!(res.len() == 2);
        assert_eq!(res["chr3"][0], (51289, 19238568));
        assert_eq!(res["chrX"][0], (346798, 23689090));
        assert_eq!(res["chrX"][1], (346798, 23689090));
    }
}
