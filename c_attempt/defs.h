const int64_t MAX_COLUMN_CHUNK_SIZE = 1024 * 1024 * 10; // 10 MB

const char* GBAM_MAGIC = "geeBAM20";

enum ColumnType {
    /// Data fields, including both constant and variable size fields.
    COLUMNTYPE_refID,
    COLUMNTYPE_pos,
    COLUMNTYPE_mapq,
    COLUMNTYPE_bin,
    COLUMNTYPE_flag,
    COLUMNTYPE_next_refID,
    COLUMNTYPE_next_pos,
    COLUMNTYPE_tlen,
    COLUMNTYPE_read_name,
    COLUMNTYPE_cigar,
    COLUMNTYPE_seq,
    COLUMNTYPE_qual,
    COLUMNTYPE_tags,
    /// Index fields (constant size), used to index where var size field recs end.
    COLUMNTYPE_index_read_name,
    COLUMNTYPE_index_cigar,
    COLUMNTYPE_index_seq,
    COLUMNTYPE_index_tags,
    COLUMNTYPE_index_qual,
};

static const char* ColumnTypeNames[] = {
    "COLUMNTYPE_refID",
    "COLUMNTYPE_pos",
    "COLUMNTYPE_mapq",
    "COLUMNTYPE_bin",
    "COLUMNTYPE_flag",
    "COLUMNTYPE_next_refID",
    "COLUMNTYPE_next_pos",
    "COLUMNTYPE_tlen",
    "COLUMNTYPE_read_name",
    "COLUMNTYPE_cigar",
    "COLUMNTYPE_seq",
    "COLUMNTYPE_qual",
    "COLUMNTYPE_tags",
    "COLUMNTYPE_index_read_name",
    "COLUMNTYPE_index_cigar",
    "COLUMNTYPE_index_seq",
    "COLUMNTYPE_index_tags",
    "COLUMNTYPE_index_qual"
};

const int64_t COLUMNTYPE_SIZE = static_cast<int64_t>(COLUMNTYPE_index_qual) + 1;

struct ColumnChunkMeta
{
    int64_t rec_num;
    int64_t file_offset;       // offset in the file where this column starts
    int64_t uncompressed_size; // size of the uncompressed data
    int64_t compressed_size;   // size of the compressed data
    ColumnChunkMeta *next;      // pointer to the next metadata chunk
    ColumnChunkMeta *prev;      // pointer to the next metadata chunk
};

struct Column {
    uint8_t *data;
    int64_t cur_ptr;
    Column* index_column;
};

struct Writer {
    int fd;
    bam_hdr_t *header;
    Column *columns;
    ColumnChunkMeta **metadatas; // Pointer to metadata for each column
    int64_t[COLUMNTYPE_SIZE] cur_chunk_rec_num; // Number of records in the current chunk for each column
};

struct Reader {
    int fd;

    int64_t rec_num;
    bam_hdr_t *header;

    Column *columns;
    ColumnChunkMeta **metadatas; // Pointer to metadata for each column

    int64_t[COLUMNTYPE_SIZE] metadatas_length; 
    int64_t[COLUMNTYPE_SIZE] currently_loaded_column_chunk; 
    int64_t[COLUMNTYPE_SIZE] loaded_since_rec_num; // The last record number loaded for each column
    int64_t[COLUMNTYPE_SIZE] loaded_up_to_rec_num; // The last record number loaded for each column
    int64_t* [COLUMNTYPE_SIZE] record_counts_per_column_chunk; 
};
