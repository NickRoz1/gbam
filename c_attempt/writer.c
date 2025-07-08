#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <htslib/sam.h>
#include <stdbool.h>
#include "utils.h"
#include "defs.h"
#include "meta.h"
#include <zlib.h>
#include <lz4.h>
#include <brotli/encode.h>
#include <json-c/json.h>
#include <zstd.h>
#include <omp.h>

#define WRITE_INDEX_COLUMN(COL_NAME)                                                                                                                 \
    write_int32_le(&columns[COL_NAME].index_column->data[columns[COL_NAME].index_column->cur_ptr], columns[COL_NAME].cur_ptr); \
    columns[COL_NAME].index_column->cur_ptr += sizeof(int32_t);

FILE *log_fp = NULL;

void init_columns(Writer *writer) {
    // Allocate memory for columns
    writer->columns = malloc(sizeof(Column) * COLUMNTYPE_SIZE); 
    writer->metadatas = malloc(sizeof(ColumnChunkMeta*) * COLUMNTYPE_SIZE); 

    if (!writer->columns) {
        fprintf(stderr, "Failed to allocate memory for columns\n");
        return;
    }

    for (int i = 0; i < COLUMNTYPE_SIZE; i++) {
        writer->columns[i].capacity = MAX_COLUMN_CHUNK_SIZE;
        writer->columns[i].data = malloc(MAX_COLUMN_CHUNK_SIZE);
        if (!writer->columns[i].data) {
            fprintf(stderr, "Failed to allocate memory for column data\n");
            return;
        }
        writer->columns[i].cur_ptr = 0;
        writer->columns[i].index_column = NULL; // Initialize index_column to NULL
        writer->metadatas[i] = NULL;
    }

    writer->columns[COLUMNTYPE_read_name].index_column = &writer->columns[COLUMNTYPE_index_read_name];
    writer->columns[COLUMNTYPE_cigar].index_column     = &writer->columns[COLUMNTYPE_index_cigar];
    writer->columns[COLUMNTYPE_seq].index_column       = &writer->columns[COLUMNTYPE_index_seq];
    writer->columns[COLUMNTYPE_qual].index_column      = &writer->columns[COLUMNTYPE_index_qual];
    writer->columns[COLUMNTYPE_tags].index_column      = &writer->columns[COLUMNTYPE_index_tags];

}

void load_codec_map(const char* json_path, Writer* writer) {
    FILE *fp = fopen(json_path, "r");
    if (!fp) {
        perror("Could not open codec map file");
        exit(1);
    }

    fseek(fp, 0, SEEK_END);
    long len = ftell(fp);
    rewind(fp);

    char *buffer = malloc(len + 1);
    fread(buffer, 1, len, fp);
    buffer[len] = '\0';
    fclose(fp);

    struct json_object *root = json_tokener_parse(buffer);
    if (!root) {
        fprintf(stderr, "Invalid JSON in codec map\n");
        free(buffer);
        exit(1);
    }

    for (int i = 0; i < COLUMNTYPE_SIZE; i++) {
        struct json_object *val;
        if (json_object_object_get_ex(root, ColumnTypeNames[i], &val)) {
            strncpy(writer->codec_map[i], json_object_get_string(val), sizeof(writer->codec_map[i]) - 1);
        } else {
            strcpy(writer->codec_map[i], "none");  // default if not specified
        }
    }

    json_object_put(root);
    free(buffer);
}

Writer* create_writer(FILE* fd, bam_hdr_t *header) {
    Writer *writer = calloc(1, sizeof(Writer));
    if (!writer) {
        fprintf(stderr, "Failed to allocate memory for writer\n");
        return NULL;
    }
    writer->fd = fd;

    log_fp = fopen("compression_log.txt", "w");
    if (!log_fp) {
        perror("Failed to open compression log file");
        exit(1);
    }
    
    // Reserve space for the header
    if (fseek(writer->fd, 1000, SEEK_SET) != 0) {
        perror("Failed to move file cursor");
        free(writer);
        return NULL;
    }

    writer->metadatas = NULL;
    writer->header = header;
    writer->columns = NULL; // Initialize columns to NULL

    load_codec_map("codec_map.json", writer);

    init_columns(writer);

    return writer;
}


int compress_buffer(const void* input_data, size_t input_size, 
                   void** output_data, size_t* output_size) {
    // Calculate maximum possible compressed size
    uLongf compressed_size = compressBound(input_size);
    
    // Allocate output buffer
    *output_data = malloc(compressed_size);
    if (!*output_data) {
        return Z_MEM_ERROR;
    }
    
    // Compress the data
    int result = compress((Bytef*)*output_data, &compressed_size,
                         (const Bytef*)input_data, input_size);
    
    if (result != Z_OK) {
        free(*output_data);
        *output_data = NULL;
        return result;
    }
    
    // Resize buffer to actual compressed size to save memory
    void* resized = realloc(*output_data, compressed_size);
    if (resized) {
        *output_data = resized;
    }
    
    *output_size = compressed_size;
    return Z_OK;
}

// int compress_buffer_brotli(const uint8_t* input, size_t input_size,
//                            uint8_t** output, size_t* output_size) {
//     size_t max_compressed_size = BrotliEncoderMaxCompressedSize(input_size);
//     *output = malloc(max_compressed_size);
//     if (*output == NULL) return -1;

//     *output_size = max_compressed_size;

//     if (!BrotliEncoderCompress(
//             5, BROTLI_DEFAULT_WINDOW, BROTLI_MODE_GENERIC,
//             input_size, input, output_size, *output)) {
//         free(*output);
//         *output = NULL;
//         return -1;
//     }

//     return 0;
// }

int compress_buffer_zstd(const void* input, size_t input_size,
                         void** output, size_t* output_size) {
    size_t bound = ZSTD_compressBound(input_size);
    *output = malloc(bound);
    if (!*output) return 1;

    size_t actual = ZSTD_compress(*output, bound, input, input_size, 3);
    if (ZSTD_isError(actual)) {
        fprintf(stderr, "Zstd compression error: %s\n", ZSTD_getErrorName(actual));
        free(*output);
        return 1;
    }
    *output_size = actual;
    return 0;
}

int compress_buffer_brotli(const uint8_t* input, size_t input_size,
                           uint8_t** output, size_t* output_size) {
    size_t max_compressed_size = BrotliEncoderMaxCompressedSize(input_size);
    uint8_t* temp_output = malloc(max_compressed_size);
    if (temp_output == NULL) return -1;

    size_t actual_size = max_compressed_size;

    BROTLI_BOOL result = BrotliEncoderCompress(
        8, BROTLI_DEFAULT_WINDOW, BROTLI_MODE_GENERIC,
        input_size, input, &actual_size, temp_output);

    if (!result) {
        free(temp_output);
        *output = NULL;
        return -1;
    }

    *output = temp_output;
    *output_size = actual_size;
    return 0;
}

int compress_data_lz4(const char* source, int source_size, 
                     char** compressed, int* compressed_size) {
    // Calculate maximum possible compressed size
    int max_compressed_size = LZ4_compressBound(source_size);
    
    // Allocate buffer for compressed data
    *compressed = malloc(max_compressed_size);
    if (!*compressed) {
        return -1;
    }
    
    // Compress the data
    *compressed_size = LZ4_compress_default(source, *compressed, 
                                           source_size, max_compressed_size);
    
    if (*compressed_size <= 0) {
        free(*compressed);
        *compressed = NULL;
        return -1;
    }
    
    // Resize to actual compressed size
    char* resized = realloc(*compressed, *compressed_size);
    if (resized) {
        *compressed = resized;
    }
    
    return 0;
}

void check_and_dump_if_full(Writer *writer, bool force_dump) {
    Column *columns = writer->columns;

    #pragma omp parallel for
    for (int i = 0; i < COLUMNTYPE_SIZE; i++) {
        if (!(force_dump || (columns[i].cur_ptr >= MAX_COLUMN_CHUNK_SIZE)))
            continue;

        ColumnChunkMeta *metadata = calloc(1, sizeof(ColumnChunkMeta));
        metadata->file_offset = ftell(writer->fd); // Still problematic in parallel context
        metadata->uncompressed_size += columns[i].cur_ptr;

        const char* codec_to_use = writer->codec_map[i];
        strncpy(metadata->codec, codec_to_use, sizeof(metadata->codec) - 1);
        metadata->codec[sizeof(metadata->codec) - 1] = '\0';

        void* compressed_data = NULL;
        size_t compressed_size = 0;
        int error = 0;

        if (strcmp(codec_to_use, "zlib") == 0) {
            error = compress_buffer(columns[i].data, columns[i].cur_ptr, &compressed_data, &compressed_size) != Z_OK;
        } else if (strcmp(codec_to_use, "lz4") == 0) {
            int int_size = 0;
            error = compress_data_lz4(columns[i].data, columns[i].cur_ptr, (char**)&compressed_data, &int_size);
            compressed_size = int_size;
        } else if (strcmp(codec_to_use, "brotli") == 0) {
            error = compress_buffer_brotli(columns[i].data, columns[i].cur_ptr, (uint8_t**)&compressed_data, &compressed_size);
        } else if (strcmp(codec_to_use, "zstd") == 0) {
            error = compress_buffer_zstd(columns[i].data, columns[i].cur_ptr, &compressed_data, &compressed_size);
        }

        if (error != 0) {
            fprintf(stderr, "Compression failed for column %d (%s)\n", i, codec_to_use);
            exit(1);
        }

        // Write compressed data and metadata back in critical section
        #pragma omp critical
        {
            ssize_t written = fwrite(compressed_data, 1, compressed_size, writer->fd);
            metadata->compressed_size += written;

            if (log_fp && columns[i].cur_ptr > 0 && strcmp(codec_to_use, "brotli") == 0) {
                fprintf(log_fp,
                        "Field: %s, Uncompressed Size: %ld, Compressed Size: %ld\n",
                        ColumnTypeNames[i],
                        metadata->uncompressed_size,
                        compressed_size);
                fflush(log_fp);
            }

            metadata->rec_num = writer->cur_chunk_rec_num[i];
            writer->cur_chunk_rec_num[i] = 0;
            columns[i].cur_ptr = 0;

            if (writer->metadatas[i] == NULL) {
                writer->metadatas[i] = metadata;
            } else {
                writer->metadatas[i]->next = metadata;
                writer->metadatas[i]->next->prev = writer->metadatas[i];
                writer->metadatas[i] = writer->metadatas[i]->next;
            }

            free(compressed_data);
        }
    }
}

int write_bam_record(Writer *writer, bam1_t *aln) {
    if (!writer || !aln) {
        fprintf(stderr, "Invalid writer or alignment record\n");
        return 1;
    }

    // Write each field to the appropriate column
    Column *columns = writer->columns;

    // Constant sized fields
    // TODO: watchout for endianess issues.
    {
        write_int32_le(&columns[COLUMNTYPE_refID].data[columns[COLUMNTYPE_refID].cur_ptr], aln->core.tid);
        columns[COLUMNTYPE_refID].cur_ptr += sizeof(int32_t);
        write_int32_le(&columns[COLUMNTYPE_pos].data[columns[COLUMNTYPE_pos].cur_ptr], aln->core.pos);
        columns[COLUMNTYPE_pos].cur_ptr += sizeof(int32_t);
        columns[COLUMNTYPE_mapq].data[columns[COLUMNTYPE_mapq].cur_ptr] = aln->core.qual;
        columns[COLUMNTYPE_mapq].cur_ptr += sizeof(uint8_t);
        write_int16_le(&columns[COLUMNTYPE_bin].data[columns[COLUMNTYPE_bin].cur_ptr], aln->core.bin);
        columns[COLUMNTYPE_bin].cur_ptr += sizeof(uint16_t);
        write_int16_le(&columns[COLUMNTYPE_flag].data[columns[COLUMNTYPE_flag].cur_ptr], aln->core.flag);
        columns[COLUMNTYPE_flag].cur_ptr += sizeof(uint16_t);
        write_int32_le(&columns[COLUMNTYPE_next_refID].data[columns[COLUMNTYPE_next_refID].cur_ptr], aln->core.mtid);
        columns[COLUMNTYPE_next_refID].cur_ptr += sizeof(int32_t);
        write_int32_le(&columns[COLUMNTYPE_next_pos].data[columns[COLUMNTYPE_next_pos].cur_ptr], aln->core.mpos);
        columns[COLUMNTYPE_next_pos].cur_ptr += sizeof(int32_t);
        write_int32_le(&columns[COLUMNTYPE_tlen].data[columns[COLUMNTYPE_tlen].cur_ptr], aln->core.isize);
        columns[COLUMNTYPE_tlen].cur_ptr += sizeof(int32_t);
    }    
    // Variable sized fields
    // TODO: handle when input cigar is bigger than MAX_COLUMN_CHUNK_SIZE (reallocate, extend buffer)
    {
        // l_extranul is not correct apparently. Strlen is slow but no other option...
        // memcpy(&columns[COLUMNTYPE_read_name].data[columns[COLUMNTYPE_read_name].cur_ptr], bam_get_qname(aln), strlen(bam_get_qname(aln)) + 1);
        // columns[COLUMNTYPE_read_name].cur_ptr += strlen(bam_get_qname(aln)) + 1; // Exclude the trailing null byte
        // int l_qname = aln->core.l_qname;
        // memcpy(&columns[COLUMNTYPE_read_name].data[columns[COLUMNTYPE_read_name].cur_ptr],
        //     bam_get_qname(aln), l_qname);
        // columns[COLUMNTYPE_read_name].cur_ptr += l_qname;
        // WRITE_INDEX_COLUMN(COLUMNTYPE_read_name)
        
        // memcpy(&columns[COLUMNTYPE_cigar].data[columns[COLUMNTYPE_cigar].cur_ptr], (char*)bam_get_cigar(aln), aln->core.n_cigar<<2);
        // columns[COLUMNTYPE_cigar].cur_ptr += aln->core.n_cigar<<2;
        // WRITE_INDEX_COLUMN(COLUMNTYPE_cigar)

        // memcpy(&columns[COLUMNTYPE_seq].data[columns[COLUMNTYPE_seq].cur_ptr], bam_get_seq(aln), (((aln)->core.l_qseq + 1)>>1));
        // columns[COLUMNTYPE_seq].cur_ptr += (((aln)->core.l_qseq + 1)>>1);
        // WRITE_INDEX_COLUMN(COLUMNTYPE_seq)

        // memcpy(&columns[COLUMNTYPE_qual].data[columns[COLUMNTYPE_qual].cur_ptr], bam_get_qual(aln), aln->core.l_qseq);
        // columns[COLUMNTYPE_qual].cur_ptr += aln->core.l_qseq;
        // WRITE_INDEX_COLUMN(COLUMNTYPE_qual)

        // memcpy(&columns[COLUMNTYPE_tags].data[columns[COLUMNTYPE_tags].cur_ptr], bam_get_aux(aln), bam_get_l_aux(aln));
        // columns[COLUMNTYPE_tags].cur_ptr += bam_get_l_aux(aln);
        // WRITE_INDEX_COLUMN(COLUMNTYPE_tags)

        int l_qname = aln->core.l_qname;
        ensure_column_capacity(&columns[COLUMNTYPE_read_name], l_qname);
        memcpy(&columns[COLUMNTYPE_read_name].data[columns[COLUMNTYPE_read_name].cur_ptr],
                bam_get_qname(aln), l_qname);
        columns[COLUMNTYPE_read_name].cur_ptr += l_qname;
        WRITE_INDEX_COLUMN(COLUMNTYPE_read_name)

        int cigar_size = aln->core.n_cigar << 2;
        ensure_column_capacity(&columns[COLUMNTYPE_cigar], cigar_size);
        memcpy(&columns[COLUMNTYPE_cigar].data[columns[COLUMNTYPE_cigar].cur_ptr],
                bam_get_cigar(aln), cigar_size);
        columns[COLUMNTYPE_cigar].cur_ptr += cigar_size;
        WRITE_INDEX_COLUMN(COLUMNTYPE_cigar)

        int seq_size = (aln->core.l_qseq + 1) >> 1;
        ensure_column_capacity(&columns[COLUMNTYPE_seq], seq_size);
        memcpy(&columns[COLUMNTYPE_seq].data[columns[COLUMNTYPE_seq].cur_ptr],
                bam_get_seq(aln), seq_size);
        columns[COLUMNTYPE_seq].cur_ptr += seq_size;
        WRITE_INDEX_COLUMN(COLUMNTYPE_seq)

        int qual_size = aln->core.l_qseq;
        ensure_column_capacity(&columns[COLUMNTYPE_qual], qual_size);
        memcpy(&columns[COLUMNTYPE_qual].data[columns[COLUMNTYPE_qual].cur_ptr],
                bam_get_qual(aln), qual_size);
        columns[COLUMNTYPE_qual].cur_ptr += qual_size;
        WRITE_INDEX_COLUMN(COLUMNTYPE_qual)


        int tag_size = bam_get_l_aux(aln);
        ensure_column_capacity(&columns[COLUMNTYPE_tags], tag_size);
        memcpy(&columns[COLUMNTYPE_tags].data[columns[COLUMNTYPE_tags].cur_ptr], bam_get_aux(aln), bam_get_l_aux(aln));
        columns[COLUMNTYPE_tags].cur_ptr += bam_get_l_aux(aln);
        WRITE_INDEX_COLUMN(COLUMNTYPE_tags)
    }


    for (int i = 0; i < COLUMNTYPE_SIZE; i++) writer->cur_chunk_rec_num[i]++;

    check_and_dump_if_full(writer, false);
    return 0;
}

void ensure_column_capacity(Column* col, int needed_size) {
    if (col->cur_ptr + needed_size <= col->capacity) {
        return; // enough space
    }

    // Double the buffer size until it's big enough
    int new_capacity = col->capacity;
    while (col->cur_ptr + needed_size > new_capacity) {
        new_capacity *= 2;
    }

    void* new_data = realloc(col->data, new_capacity);
    if (!new_data) {
        fprintf(stderr, "Failed to realloc buffer to %d bytes\n", new_capacity);
        exit(1);
    }

    col->data = new_data;
    col->capacity = new_capacity;
}


void close_writer(Writer *writer) {
    if (!writer) return;

    // Write any remaining data in columns
    check_and_dump_if_full(writer, true);
    const int64_t cur_file_offset = ftell(writer->fd);
    write_meta(writer->fd, writer->metadatas, COLUMNTYPE_SIZE);
    fprintf(writer->fd, "\0");
    int64_t meta_size = ftell(writer->fd) - cur_file_offset;

    // Write SAM header after the metadata
    int32_t header_len = sam_hdr_length(writer->header);
    fwrite(&header_len, sizeof(int32_t), 1, writer->fd);
    fwrite(sam_hdr_str(writer->header), 1, header_len, writer->fd);

    // Seek beginning of the file to write the header
    fseek(writer->fd, 0, SEEK_SET);
    write_header(writer->fd, cur_file_offset, meta_size, header_len);

    // Free allocated memory for columns and metadata
    for (int i = 0; i < COLUMNTYPE_SIZE; i++) {
        free(writer->columns[i].data);

        // Go to head of metadata list
        ColumnChunkMeta *meta = writer->metadatas[i];
        while (meta && meta->prev) {
            meta = meta->prev;
        }

        // Free entire metadata chain
        while (meta) {
            ColumnChunkMeta *tmp = meta;
            meta = meta->next;
            free(tmp);
        }
    }

    free(writer->columns);
    free(writer->metadatas);
    free(writer);

    if (log_fp) {
        fclose(log_fp);
        log_fp = NULL;
    }
}
