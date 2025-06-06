#pragma once

#include <stdint.h>
#include <htslib/sam.h>  
  
static const int64_t MAX_COLUMN_CHUNK_SIZE = 1024 * 1024 * 10; // 10 MB

static const char* GBAM_MAGIC = "geeBAM20";

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
    COLUMNTYPE_SIZE // Total number of column types
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

typedef struct ColumnChunkMeta
{
    int64_t rec_num;
    int64_t file_offset;       // offset in the file where this column starts
    int64_t uncompressed_size; // size of the uncompressed data
    int64_t compressed_size;   // size of the compressed data
    struct ColumnChunkMeta *next;      // pointer to the next metadata chunk
    struct ColumnChunkMeta *prev;      // pointer to the next metadata chunk
} ColumnChunkMeta;

typedef struct Column
{
    uint8_t *data;
    int64_t cur_ptr;
    struct Column* index_column;
} Column;

typedef struct  {
    FILE* fd;
    bam_hdr_t *header;
    Column *columns;
    ColumnChunkMeta **metadatas; // Pointer to metadata for each column
    int64_t cur_chunk_rec_num[COLUMNTYPE_SIZE] ; // Number of records in the current chunk for each column
} Writer;

typedef struct {
    int fd;

    int64_t rec_num;
    bam_hdr_t *header;

    Column *columns;
    ColumnChunkMeta **metadatas; // Pointer to metadata for each column

    int64_t metadatas_lengths[COLUMNTYPE_SIZE]; 
    int64_t currently_loaded_column_chunk[COLUMNTYPE_SIZE]; 
    int64_t loaded_since_rec_num[COLUMNTYPE_SIZE]; // The last record number loaded for each column
    int64_t loaded_up_to_rec_num[COLUMNTYPE_SIZE]; // The last record number loaded for each column
    int64_t*  record_counts_per_column_chunk[COLUMNTYPE_SIZE]; 
} Reader;

struct bam1_t;
struct bam_hdr_t ;

Writer *create_writer(FILE* fd, bam_hdr_t  *header);
int write_bam_record(Writer *writer, bam1_t *aln);
void close_writer(Writer *writer);

Reader *make_reader(FILE *fp);
void read_record(Reader *reader, int64_t rec_num, bam1_t *aln);