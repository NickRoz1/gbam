use super::u8_size;
use std::io::Write;
/// Type of compression to utilize in writer.
/// Serialized into two bytes - type of compression and its strength.
#[allow(missing_docs)]
#[derive(PartialEq, Debug, Clone)]
pub enum Compression {
    ZSTD(u8),
    FLATE2(u8),
}

impl Default for Compression {
    fn default() -> Self {
        Compression::ZSTD(255)
    }
}

pub const COMPRESSION_ENUM_SIZE: usize = u8_size * 2;

impl From<&[u8]> for Compression {
    fn from(mut bytes: &[u8]) -> Self {
        assert!(
            bytes.len() >= COMPRESSION_ENUM_SIZE,
            "Not enough bytes to form column chunk meta struct.",
        );
        let compr = bytes[0];
        match compr {
            0 => {
                // ZSTD
                let strength = bytes[1];
                if strength < 1 || strength > 21 {
                    panic!("Metadata is damaged. ZSTD compression strength should be in limits [1, 21].");
                }
                return Compression::ZSTD(strength);
            }
            1 => {
                let strength = bytes[1];
                if strength < 0 || strength > 9 {
                    panic!("Metadata is damaged. FLATE2 compression strength should be in limits [0, 9].");
                }
                return Compression::FLATE2(strength);
            }
            _ => panic!("This GBAM version doesn't support this compression algorithm."),
        }
    }
}

impl Into<Vec<u8>> for Compression {
    fn into(self) -> Vec<u8> {
        let mut res = Vec::<u8>::new();
        self.dump_as_bytes(&mut res).unwrap();
        res
    }
}

impl Compression {
    pub fn dump_as_bytes<T: Write>(&self, writer: &mut T) -> std::io::Result<usize> {
        match self {
            Compression::ZSTD(strength) => writer.write(&[0, strength.to_owned()]),
            Compression::FLATE2(strength) => writer.write(&[1, strength.to_owned()]),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_from() -> () {
        let zstd_correct = &[0, 21];
        let flate2_correct = &[1, 2];

        let mut compr = Compression::from(&zstd_correct[..]);
        assert_eq!(compr, Compression::ZSTD(21));
        compr = Compression::from(&flate2_correct[..]);
        assert_eq!(compr, Compression::FLATE2(2));
    }

    #[test]
    fn test_into() -> () {
        let zstd_correct = &[0, 21];
        let flate2_correct = &[1, 2];

        let zstd = Compression::ZSTD(21);
        let flate2 = Compression::FLATE2(2);

        let bytes_zstd: Vec<u8> = zstd.into();
        let bytes_flate2: Vec<u8> = flate2.into();

        assert_eq!(&bytes_zstd[..], zstd_correct);
        assert_eq!(&bytes_flate2[..], flate2_correct);
    }

    #[test]
    #[should_panic(expected = "This GBAM version doesn't support this compression algorithm.")]
    fn test_unknown_compression_panic() {
        let incorrect = &[200, 1];
        Compression::from(&incorrect[..]);
    }

    #[test]
    #[should_panic(expected = "Metadata is damaged. ZSTD compression strength should be in limits")]
    fn test_wrong_data() {
        let zstd_incorrect = &[0, 22];
        Compression::from(&zstd_incorrect[..]);
    }
}
