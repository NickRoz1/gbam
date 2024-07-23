#include "../api.h"
#include <dlfcn.h>
#include <assert.h>


int main(){
    void* reader = get_gbam_reader("/home/nickr/projects/gbam/test_data/one_record.gbam");
    sam_hdr_t* hdr = get_header(reader);
    // printf("%s", sam_hdr_str(hdr));
    char* my_str = (char*) calloc(10000*sizeof(char), 0);
    kstring_t s = {10000,10000,my_str};
    bam1_t rec = get_bam_record(reader, 0);
    sam_format1(hdr, &rec, &s);
    printf("%d", rec.core.n_cigar);
    printf("%s", s.s);
    // free_header(hdr);
}