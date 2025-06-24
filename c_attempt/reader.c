#include "defs.h"
#include <json-c/json.h>
#include "utils.h"
#include "assert.h"
#include <htslib/hts.h>
#include <stdbool.h>
#include <zlib.h>
#include <lz4.h>


#define ADJUSTED_OFFSET(COLUMNTYPE) \
    (rec_num-(reader->loaded_since_rec_num[COLUMNTYPE]))

#define bam_reg2bin(beg,end) hts_reg2bin((beg),(end),14,5)

void parse_meta_from_json_string(const char *json_str, Reader *reader);


static void bam_cigar2rqlens(int n_cigar, const uint32_t *cigar,
                             hts_pos_t *rlen, hts_pos_t *qlen)
{
    int k;
    *rlen = *qlen = 0;
    for (k = 0; k < n_cigar; ++k) {
        int type = bam_cigar_type(bam_cigar_op(cigar[k]));
        int len = bam_cigar_oplen(cigar[k]);
        if (type & 1) *qlen += len;
        if (type & 2) *rlen += len;
    }
}


Reader* make_reader(char* file){
    if (!file) {
        fprintf(stderr, "File pointer is NULL\n");
        return NULL;
    }

    struct json_object *root = json_tokener_parse(file);

    int64_t seekpos = -1;
    int64_t meta_size = -1;

    struct json_object *seekpos_obj;
    if (json_object_object_get_ex(root, "seekpos", &seekpos_obj)) {
        seekpos = json_object_get_int64(seekpos_obj);
    }
    struct json_object *meta_size_obj;
    if (json_object_object_get_ex(root, "meta_size", &meta_size_obj)) {
        meta_size = json_object_get_int64(meta_size_obj);
    }

    json_object_put(root); 


    assert(seekpos >= 0);
    assert(meta_size >= 0);

    Reader* reader = (Reader*)calloc(1, sizeof(Reader));

    parse_meta_from_json_string(file+seekpos, reader);

    reader->mmaped_file = file;
    reader->header = sam_hdr_parse(strlen(file+seekpos+meta_size), file+seekpos+meta_size);
    reader->rec_num = 0;
    reader->columns = (Column*)calloc(COLUMNTYPE_SIZE, sizeof(Column));

    for(int i = 0; i < COLUMNTYPE_SIZE; i++) {
        ColumnChunkMeta *metas = reader->metadatas[i];
        reader->record_counts_per_column_chunk[i] = (int64_t*)calloc(reader->metadatas_lengths[i], sizeof(int64_t));
        int64_t accum = 0;
        for (int j = 0; j < reader->metadatas_lengths[i]; j++)
        {
            accum += metas[j].rec_num;
            reader->record_counts_per_column_chunk[i][j] = accum;
        }
        reader->rec_num = accum;
    }

    for(int i = 0; i < COLUMNTYPE_SIZE; i++) {
        reader->currently_loaded_column_chunk[i] = -1;
        reader->loaded_up_to_rec_num[i] = -1;
        reader->loaded_since_rec_num[i] = -1;
    }

    return reader;
}

int decompress_buffer(const void* compressed_data, size_t compressed_size,
                     void** output_data, uLongf decompressed_size) {
    int result = uncompress((Bytef*)*output_data, &decompressed_size,
                           (const Bytef*)compressed_data, compressed_size);
    
    if (result != Z_OK) {
        free(*output_data);
        *output_data = NULL;
        return result;
    }
    
    return Z_OK;
}


void fetch_field(Reader* reader, int64_t rec_num, int64_t COLUMNTYPE){
    if(reader->loaded_up_to_rec_num[COLUMNTYPE] >= rec_num && 
       reader->loaded_since_rec_num[COLUMNTYPE] <= rec_num) {
        return;
    }

    int l = -1;
    int r = reader->metadatas_lengths[COLUMNTYPE]-1;
    
    while((r-l)>1){
        int m = (l+r)/2;
        if (rec_num <= reader->record_counts_per_column_chunk[COLUMNTYPE][m]) {
            l = m;
        } else {
            r = m;
        }
    }
    if(r != reader->currently_loaded_column_chunk[COLUMNTYPE]) {
        reader->currently_loaded_column_chunk[COLUMNTYPE] = r;
        ColumnChunkMeta *meta = &reader->metadatas[COLUMNTYPE][r];
        if(reader->m_chunk_memory[COLUMNTYPE] < meta->uncompressed_size){
            if(reader->columns[COLUMNTYPE].data) {
                free(reader->columns[COLUMNTYPE].data);
            }
            reader->columns[COLUMNTYPE].data = (uint8_t*)malloc(meta->uncompressed_size);
            reader->m_chunk_memory[COLUMNTYPE] = meta->uncompressed_size;
        }
        if (!reader->columns[COLUMNTYPE].data) {
            fprintf(stderr, "Failed to allocate memory for column data\n");
            return;
        }
        void* read_buffer = reader->mmaped_file+meta->file_offset;

        if(strcmp(meta->codec, "zlib") == 0){
            int res = decompress_buffer(read_buffer, meta->compressed_size, 
                    (void**)&reader->columns[COLUMNTYPE].data, 
                    meta->uncompressed_size);
            if(res != Z_OK){
                fprintf(stderr, "Failed to decompress column data: %d\n", res);
                exit(1);
            }   
        }
        else if(strcmp(meta->codec, "lz4") == 0){
            int decompressed_size = LZ4_decompress_safe(read_buffer,
                                    reader->columns[COLUMNTYPE].data,
                                    meta->compressed_size, 
                                    meta->uncompressed_size);
            if (decompressed_size < 0) {
                printf("LZ4 decompression failed with error code: %d\n", decompressed_size);
                return;
            }
        }
        else{
            // No compression, just copy the data
            memcpy(reader->columns[COLUMNTYPE].data, read_buffer, meta->uncompressed_size);
        }

        if(r == 0){
            reader->loaded_since_rec_num[COLUMNTYPE] = 0;
        }
        else{
            reader->loaded_since_rec_num[COLUMNTYPE] = reader->record_counts_per_column_chunk[COLUMNTYPE][r-1];
        }
        reader->loaded_up_to_rec_num[COLUMNTYPE] = reader->record_counts_per_column_chunk[COLUMNTYPE][r]-1;
    }
    else{
        // Impossible to reach here, but just in case.
        exit(1);
    }
}

void preload_chunks(Reader* reader, int64_t rec_num) {
    for(int i = 0; i < COLUMNTYPE_SIZE; i++) {
        fetch_field(reader, rec_num, i);
    }
}//

void read_record(Reader* reader, int64_t rec_num, bam1_t* aln) {
    preload_chunks(reader, rec_num);

    int32_t refID = read_int32_le(&reader->columns[COLUMNTYPE_refID].data[ADJUSTED_OFFSET(COLUMNTYPE_refID) * sizeof(int32_t)]);
    int32_t pos = read_int32_le(&reader->columns[COLUMNTYPE_pos].data[ADJUSTED_OFFSET(COLUMNTYPE_pos) * sizeof(int32_t)]);
    int8_t mapq = reader->columns[COLUMNTYPE_mapq].data[ADJUSTED_OFFSET(COLUMNTYPE_mapq)];
    int16_t bin = read_int16_le(&reader->columns[COLUMNTYPE_bin].data[ADJUSTED_OFFSET(COLUMNTYPE_bin) * sizeof(int16_t)] );
    int16_t flag = read_int16_le(&reader->columns[COLUMNTYPE_flag].data[ADJUSTED_OFFSET(COLUMNTYPE_flag) * sizeof(int16_t)] );
    int32_t next_pos = read_int32_le(&reader->columns[COLUMNTYPE_next_pos].data[ADJUSTED_OFFSET(COLUMNTYPE_next_pos) * sizeof(int32_t)]);
    int32_t next_refID = read_int32_le(&reader->columns[COLUMNTYPE_next_refID].data[ADJUSTED_OFFSET(COLUMNTYPE_next_refID) * sizeof(int32_t)] );
    int32_t tlen = read_int32_le(&reader->columns[COLUMNTYPE_tlen].data[ADJUSTED_OFFSET(COLUMNTYPE_tlen) * sizeof(int32_t)]);
    
    int64_t read_name_end = read_int32_le(&reader->columns[COLUMNTYPE_index_read_name].data[ADJUSTED_OFFSET(COLUMNTYPE_index_read_name) * sizeof(int32_t)]);
    int64_t read_name_beg = 0;
    if (rec_num != 0)
    {
        fetch_field(reader, rec_num - 1, COLUMNTYPE_index_read_name);
        read_name_beg = read_int32_le(&reader->columns[COLUMNTYPE_index_read_name].data[(rec_num-1-(reader->loaded_since_rec_num[COLUMNTYPE_index_read_name]))
 * sizeof(int32_t)]);
    }
    int64_t read_cigar_end = read_int32_le(&reader->columns[COLUMNTYPE_index_cigar].data[ADJUSTED_OFFSET(COLUMNTYPE_index_cigar) * sizeof(int32_t)]);
    int64_t read_cigar_beg = 0;
    if (rec_num != 0)
    {
        fetch_field(reader, rec_num - 1, COLUMNTYPE_index_cigar);
        read_cigar_beg = read_int32_le(&reader->columns[COLUMNTYPE_index_cigar].data[(rec_num-1-(reader->loaded_since_rec_num[COLUMNTYPE_index_cigar]))
 * sizeof(int32_t)]);
    }
    int64_t read_seq_end = read_int32_le(&reader->columns[COLUMNTYPE_index_seq].data[ADJUSTED_OFFSET(COLUMNTYPE_index_seq) * sizeof(int32_t)]);
    int64_t read_seq_beg = 0;
    if (rec_num != 0)
    {
        fetch_field(reader, rec_num - 1, COLUMNTYPE_index_seq);
        read_seq_beg = read_int32_le(&reader->columns[COLUMNTYPE_index_seq].data[(rec_num-1-(reader->loaded_since_rec_num[COLUMNTYPE_index_seq]))
 * sizeof(int32_t)]);
    }
    int64_t read_qual_end = read_int32_le(&reader->columns[COLUMNTYPE_index_qual].data[ADJUSTED_OFFSET(COLUMNTYPE_index_qual) * sizeof(int32_t)]);
    int64_t read_qual_beg = 0;
    if (rec_num != 0)
    {
        fetch_field(reader, rec_num - 1, COLUMNTYPE_index_qual);
        read_qual_beg = read_int32_le(&reader->columns[COLUMNTYPE_index_qual].data[(rec_num-1-(reader->loaded_since_rec_num[COLUMNTYPE_index_qual]))
 * sizeof(int32_t)]);
    }
    int64_t read_tags_end = read_int32_le(&reader->columns[COLUMNTYPE_index_tags].data[ADJUSTED_OFFSET(COLUMNTYPE_index_tags) * sizeof(int32_t)]);
    int64_t read_tags_beg = 0;
    if (rec_num != 0)
    {
        fetch_field(reader, rec_num - 1, COLUMNTYPE_index_tags);
        read_tags_beg = read_int32_le(&reader->columns[COLUMNTYPE_index_tags].data[(rec_num-1-(reader->loaded_since_rec_num[COLUMNTYPE_index_tags]))
 * sizeof(int32_t)]);
    }

    int64_t l_qname = read_name_end - read_name_beg;
    const int64_t l_cigar = read_cigar_end - read_cigar_beg;
    const int64_t l_qual = read_qual_end - read_qual_beg;
    const int64_t l_seq = read_seq_end - read_seq_beg;
    assert(l_seq == ((l_qual+1) >> 1)); // l_seq is always half of l_qual + 1, as per BAM format
    const int64_t l_tags = read_tags_end - read_tags_beg;

    size_t qname_nuls;

    // use a default qname "*" if none is provided
    bool default_name = false;
    if (l_qname == 0)
    {
        default_name = true;
        l_qname = 1;
        qname_nuls = 4 - l_qname % 4;
    }
    else{
        qname_nuls = 4 - l_qname % 4;
    }

    uint64_t bytes_we_need = l_qname + qname_nuls + l_cigar + l_seq + l_qual + l_tags;

    if(aln->m_data < bytes_we_need) {
        if (aln->data) {
            free(aln->data);
        }
        aln->data = (uint8_t*)malloc(bytes_we_need);
        aln->m_data = bytes_we_need;
    }
    aln->l_data = bytes_we_need;

    if(default_name){
        aln->data[0] = '*';
        aln->data[1] = '\0';
        aln->data[2] = '\0';
        aln->data[3] = '\0';
    }
    else{
        memcpy(aln->data, 
           &reader->columns[COLUMNTYPE_read_name].data[read_name_beg], l_qname);
        for (int i = 0; i < qname_nuls; i++) 
            aln->data[l_qname + i] = '\0'; // Fill with null bytes
    }
    
    memcpy(&aln->data[l_qname + qname_nuls], 
           &reader->columns[COLUMNTYPE_cigar].data[read_cigar_beg], l_cigar);
    memcpy(&aln->data[l_qname + qname_nuls+l_cigar], 
           &reader->columns[COLUMNTYPE_seq].data[read_seq_beg], l_seq);
    memcpy(&aln->data[l_qname + qname_nuls+l_cigar+l_seq], 
           &reader->columns[COLUMNTYPE_qual].data[read_qual_beg], l_qual);
    memcpy(&aln->data[l_qname + qname_nuls+l_cigar+l_seq+l_qual], 
           &reader->columns[COLUMNTYPE_tags].data[read_tags_beg], l_tags);

    hts_pos_t rlen = 1, qlen = 0;
    if (!(flag & BAM_FUNMAP)) {
        bam_cigar2rqlens((int)(l_cigar >> 2), &reader->columns[COLUMNTYPE_cigar].data[read_cigar_beg], &rlen, &qlen);
    }
    if (rlen == 0) {
        rlen = 1;
    }

    aln->core.pos = pos;
    aln->core.tid = refID;
    aln->core.bin = bam_reg2bin(pos, pos + rlen);
    aln->core.qual = mapq;
    aln->core.l_extranul = (uint8_t)(qname_nuls - 1);
    aln->core.flag = flag;
    aln->core.l_qname = (uint16_t)(l_qname + qname_nuls);
    aln->core.n_cigar = (uint32_t)(l_cigar >> 2); // l_cigar is in bytes, n_cigar is in 32-bit words
    aln->core.l_qseq = (int32_t)l_qual;
    aln->core.mtid = next_refID;
    aln->core.mpos = next_pos;
    aln->core.isize = tlen;
}

void close_reader(Reader *reader) {
    bam_hdr_destroy(reader->header);

    for (int i = 0; i < COLUMNTYPE_SIZE; i++) {
        free(reader->columns[i].data);
        free(reader->record_counts_per_column_chunk[i]);
        free(reader->metadatas[i]);
    }
    free(reader->columns);
    free(reader->metadatas); 
    free(reader);
    reader = NULL;
}
