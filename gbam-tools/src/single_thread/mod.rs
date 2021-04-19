pub mod meta;
pub mod reader;
pub mod writer;

const SIZE_LIMIT: usize = 16777216;
static GBAM_MAGIC: &[u8] = b"geeBAM10";
