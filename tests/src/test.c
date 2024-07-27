#include "../inc/api.h"
#include <dlfcn.h>
#include <assert.h>
#include <unistd.h>


int main(int argc, char **argv){
    
    void* reader = get_gbam_reader(argv[1]);
    sam_hdr_t* hdr = get_header(reader);
    
    for(int i = 0; i < get_records_num(reader); i++){
        kstring_t s = {0,0,NULL};

        bam1_t rec = get_bam_record(reader, i);
        int res = sam_format1(hdr, &rec, &s);
        
        printf("%s\n", s.s);

        free(s.s);
        free_bam_record(rec);
    }
    
    free_header(hdr);
}