#include "sam.h"
// api.h
#ifndef API_H
#define API_H

#ifdef __cplusplus
extern "C" {
#endif

void* get_gbam_reader(const char* gbam_path);
bam1_t get_bam_record(void* reader, uint64_u rec_num);
sam_hdr_t* get_header(void* reader);

void free_bam_record(bam1_t rec);
void free_gbam_reader(void* reader);
void free_header(sam_hdr_t* hdr);

#ifdef __cplusplus
}
#endif

#endif // API_H