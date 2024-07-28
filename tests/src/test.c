#include "../inc/api.h"
#include <dlfcn.h>
#include <assert.h>
#include <unistd.h>

int main(int argc, char **argv){
    
    void* reader = get_gbam_reader(argv[1]);
    sam_hdr_t* hdr = get_header(reader);

    bam1_t buf = get_empty_bam1_t();
    for(int i = 0; i < get_records_num(reader); i++){
        kstring_t s = {0,0,NULL};
        refill_bam_record(reader, i, &buf);
        int res = sam_format1(hdr, &buf, &s);
  
        printf("%s\n", s.s);

        free(s.s);
    }
    free_bam_record(buf);
    free_header(hdr);
}