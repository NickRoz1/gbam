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

#define WRITE_INDEX_COLUMN(COL_NAME)                                                                                                                 \
    write_int32_le(&columns[COL_NAME].index_column->data[columns[COL_NAME].index_column->cur_ptr], columns[COL_NAME].cur_ptr); \
    columns[COL_NAME].index_column->cur_ptr += sizeof(int32_t);

void init_columns(Writer *writer) {
    // Allocate memory for columns
    writer->columns = malloc(sizeof(Column) * COLUMNTYPE_SIZE); 
    writer->metadatas = malloc(sizeof(ColumnChunkMeta*) * COLUMNTYPE_SIZE); 

    if (!writer->columns) {
        fprintf(stderr, "Failed to allocate memory for columns\n");
        return;
    }

    for (int i = 0; i < COLUMNTYPE_SIZE; i++) {
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

Writer* create_writer(FILE* fd, bam_hdr_t *header) {
    Writer *writer = calloc(1, sizeof(Writer));
    if (!writer) {
        fprintf(stderr, "Failed to allocate memory for writer\n");
        return NULL;
    }
    writer->fd = fd;
    
    // Reserve space for the header
    if (fseek(writer->fd, 1000, SEEK_SET) != 0) {
        perror("Failed to move file cursor");
        free(writer);
        return NULL;
    }

    writer->metadatas = NULL;
    writer->header = header;
    writer->columns = NULL; // Initialize columns to NULL

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

int compress_buffer_brotli(const uint8_t* input, size_t input_size,
                           uint8_t** output, size_t* output_size) {
    size_t max_compressed_size = BrotliEncoderMaxCompressedSize(input_size);
    uint8_t* temp_output = malloc(max_compressed_size);
    if (temp_output == NULL) return -1;

    size_t actual_size = max_compressed_size;

    BROTLI_BOOL result = BrotliEncoderCompress(
        5, BROTLI_DEFAULT_WINDOW, BROTLI_MODE_GENERIC,
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
    for (int i = 0; i < COLUMNTYPE_SIZE; i++) {

        if (force_dump || (columns[i].cur_ptr >= MAX_COLUMN_CHUNK_SIZE))
        {

            ColumnChunkMeta *metadata = calloc(1, sizeof(ColumnChunkMeta));
            metadata->file_offset = ftell(writer->fd);
            metadata->uncompressed_size += columns[i].cur_ptr;

            char* codec_to_use = "brotli";


            memcpy(metadata->codec, codec_to_use, strlen(codec_to_use) + 1);

            ssize_t written = 0;

            if(strcmp(metadata->codec, "zlib") == 0){
                void* compressed_data = NULL;
                size_t compressed_size = 0;
                // Compress the data
                if(compress_buffer(columns[i].data, columns[i].cur_ptr, 
                    &compressed_data, &compressed_size) != Z_OK) {
                    fprintf(stderr, "Failed to compress column data\n");
                    exit(1);
                }
            
                // Write the column data to the file
                written = fwrite(compressed_data, 1, compressed_size, writer->fd);
                free(compressed_data);
            }
            else if(strcmp(metadata->codec, "lz4") == 0){
                void* compressed_data = NULL;
                int compressed_size = 0;
                
                // Compress the data
                if(compress_data_lz4(columns[i].data, columns[i].cur_ptr, 
                    &compressed_data, &compressed_size) != 0) {
                    fprintf(stderr, "Failed to compress column data\n");
                    exit(1);
                }
            
                // Write the column data to the file
                written = fwrite(compressed_data, 1, compressed_size, writer->fd);
                free(compressed_data);
            }
            else if(strcmp(metadata->codec, "brotli") == 0){
                void* compressed_data = NULL;
                size_t compressed_size = 0;

                if (compress_buffer_brotli(columns[i].data, columns[i].cur_ptr,
                                        (uint8_t**)&compressed_data, &compressed_size) != 0) {
                    fprintf(stderr, "Failed to compress column data with Brotli\n");
                    exit(1);
                }

                written = fwrite(compressed_data, 1, compressed_size, writer->fd);
                free(compressed_data);
            }

            else{
                written = fwrite(columns[i].data, 1, columns[i].cur_ptr, writer->fd);
            }

            metadata->compressed_size += written;
            metadata->rec_num = writer->cur_chunk_rec_num[i];
            writer->cur_chunk_rec_num[i] = 0; // Reset the record count for this column
            columns[i].cur_ptr = 0; // Reset the column data
            
            // Update metadata for this column
            if (writer->metadatas[i] == NULL) {
                writer->metadatas[i] = metadata; // First metadata for this column
            }
            else{
                writer->metadatas[i]->next = metadata;
                writer->metadatas[i]->next->prev = writer->metadatas[i];
                writer->metadatas[i] = writer->metadatas[i]->next;
            }
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
        memcpy(&columns[COLUMNTYPE_read_name].data[columns[COLUMNTYPE_read_name].cur_ptr], bam_get_qname(aln), strlen(bam_get_qname(aln)));
        columns[COLUMNTYPE_read_name].cur_ptr += strlen(bam_get_qname(aln)); // Exclude the trailing null byte
        WRITE_INDEX_COLUMN(COLUMNTYPE_read_name)
        
        memcpy(&columns[COLUMNTYPE_cigar].data[columns[COLUMNTYPE_cigar].cur_ptr], (char*)bam_get_cigar(aln), aln->core.n_cigar<<2);
        columns[COLUMNTYPE_cigar].cur_ptr += aln->core.n_cigar<<2;
        WRITE_INDEX_COLUMN(COLUMNTYPE_cigar)

        memcpy(&columns[COLUMNTYPE_seq].data[columns[COLUMNTYPE_seq].cur_ptr], bam_get_seq(aln), (((aln)->core.l_qseq + 1)>>1));
        columns[COLUMNTYPE_seq].cur_ptr += (((aln)->core.l_qseq + 1)>>1);
        WRITE_INDEX_COLUMN(COLUMNTYPE_seq)

        memcpy(&columns[COLUMNTYPE_qual].data[columns[COLUMNTYPE_qual].cur_ptr], bam_get_qual(aln), aln->core.l_qseq);
        columns[COLUMNTYPE_qual].cur_ptr += aln->core.l_qseq;
        WRITE_INDEX_COLUMN(COLUMNTYPE_qual)

        memcpy(&columns[COLUMNTYPE_tags].data[columns[COLUMNTYPE_tags].cur_ptr], bam_get_aux(aln), bam_get_l_aux(aln));
        columns[COLUMNTYPE_tags].cur_ptr += bam_get_l_aux(aln);
        WRITE_INDEX_COLUMN(COLUMNTYPE_tags)
    }


    for (int i = 0; i < COLUMNTYPE_SIZE; i++) writer->cur_chunk_rec_num[i]++;

    check_and_dump_if_full(writer, false);
    return 0;
}

void close_writer(Writer *writer) {
    if (!writer) return;

    // Write any remaining data in columns
    check_and_dump_if_full(writer, true);
    const int64_t cur_file_offset = ftell(writer->fd);
    write_meta(writer->fd, writer->metadatas, COLUMNTYPE_SIZE);
    fprintf(writer->fd, "\0");
    int64_t meta_size = ftell(writer->fd)-cur_file_offset;
    // Write SAM header after the metadata
    int32_t header_len = sam_hdr_length(writer->header);  // Get header length
    fwrite(&header_len, sizeof(int32_t), 1, writer->fd);  // Write length prefix
    fwrite(sam_hdr_str(writer->header), 1, header_len, writer->fd);  // Write actual header
    // TODO: Save crc32 also and everything we had in GBAM..
    // Seek beginning of the file to write the header
    fseek(writer->fd, 0, SEEK_SET);
    write_header(writer->fd, cur_file_offset, meta_size, header_len);

    // Free allocated memory for columns and metadata
    for (int i = 0; i < COLUMNTYPE_SIZE; i++) {
        free(writer->columns[i].data);
        if (writer->metadatas[i]) {
            free(writer->metadatas[i]);
        }
    }
    
    free(writer->columns);
    free(writer->metadatas);
    
    // Free the writer itself
    free(writer);
}