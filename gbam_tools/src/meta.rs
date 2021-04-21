use super::{u64_size, Compression, COMPRESSION_ENUM_SIZE, FIELDS_NUM};
use byteorder::{LittleEndian, ReadBytesExt, WriteBytesExt};
use std::io::Write;
/// File starts with magic number. The BAM table is then divided into rowgroups,
/// each rowgroup is written into the file column by column (column chunk inside
/// the rowgroup). Each column chunk compressed separately. Information about
/// column chunk offset and length is stored in the file metadata.
#[derive(Clone, Debug)]
pub struct RowGroupMeta {
    pub offset: u64,
    pub cols: Vec<ColChunkMeta>,
}

const ROW_GROUP_META_SIZE: usize = u64_size + FIELDS_NUM * COL_CHUNK_META_SIZE;

impl From<&[u8]> for RowGroupMeta {
    fn from(mut bytes: &[u8]) -> Self {
        assert!(
            bytes.len() >= ROW_GROUP_META_SIZE,
            "Not enough bytes to form row group meta struct.",
        );
        RowGroupMeta {
            offset: bytes
                .read_u64::<LittleEndian>()
                .expect("Rowgroup metadata is damaged."),
            cols: (0..FIELDS_NUM).map(|_| ColChunkMeta::from(bytes)).collect(),
        }
    }
}

impl Into<Vec<u8>> for RowGroupMeta {
    fn into(self) -> Vec<u8> {
        let mut res = Vec::<u8>::new();
        self.dump_as_bytes(&mut res);
        res
    }
}

impl RowGroupMeta {
    pub fn dump_as_bytes<T: Write>(&self, writer: &mut T) {
        writer.write_u64::<LittleEndian>(self.offset).unwrap();
        for col in &self.cols {
            col.dump_as_bytes(writer);
        }
    }
}

const COL_CHUNK_META_SIZE: usize = u64_size * 4 + COMPRESSION_ENUM_SIZE;
#[derive(Clone, Debug, PartialEq, Default)]
pub struct ColChunkMeta {
    pub offset: u64,
    pub compr_size: u64,
    pub uncompr_size: u64,
    pub val_num: u64,
    pub compressor: Compression,
}

impl From<&[u8]> for ColChunkMeta {
    fn from(mut bytes: &[u8]) -> Self {
        assert!(
            bytes.len() >= COL_CHUNK_META_SIZE,
            "Not enough bytes to form column chunk meta struct.",
        );
        ColChunkMeta {
            offset: bytes
                .read_u64::<LittleEndian>()
                .expect("Column metadata is damaged."),
            compr_size: bytes
                .read_u64::<LittleEndian>()
                .expect("Column metadata is damaged."),
            uncompr_size: bytes
                .read_u64::<LittleEndian>()
                .expect("Column metadata is damaged."),
            val_num: bytes
                .read_u64::<LittleEndian>()
                .expect("Column metadata is damaged."),
            compressor: Compression::from(bytes),
        }
    }
}

impl ColChunkMeta {
    pub fn dump_as_bytes<T: Write>(&self, writer: &mut T) {
        writer.write_u64::<LittleEndian>(self.offset).unwrap();
        writer.write_u64::<LittleEndian>(self.compr_size).unwrap();
        writer.write_u64::<LittleEndian>(self.uncompr_size).unwrap();
        writer.write_u64::<LittleEndian>(self.val_num).unwrap();
        self.compressor.dump_as_bytes(writer).unwrap();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_row_and_col_meta() -> () {
        let mut cols_info = Vec::<ColChunkMeta>::new();
        for _ in 0..FIELDS_NUM {
            cols_info.push(ColChunkMeta {
                offset: 1,
                compr_size: 2,
                uncompr_size: 3,
                val_num: 4,
                compressor: Compression::ZSTD(2),
            });
        }

        let original = RowGroupMeta {
            offset: 5,
            cols: cols_info,
        };

        let buf: Vec<u8> = original.clone().into();

        let restored = RowGroupMeta::from(&buf[..]);

        assert_eq!(restored.offset, original.offset);
        assert_eq!(&restored.cols[..], &original.cols[..]);
    }
}
