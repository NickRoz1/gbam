#include <htslib/sam.h>
#include <stdlib.h>
#include "defs.h"

int main() {
    // Initialize input from stdin
    // samFile *in = sam_open("-", "r");
    // if (!in) {
    //     fprintf(stderr, "Failed to open stdin\n");
    //     return 1;
    // }

    // // Read header
    // bam_hdr_t *header = sam_hdr_read(in);
    // if (!header) {
    //     fprintf(stderr, "Failed to read header\n");
    //     sam_close(in);
    //     return 1;
    // }

    // // Initialize output to stdout (SAM format)
    // samFile *out = sam_open("-", "w");
    // if (!out) {
    //     fprintf(stderr, "Failed to open stdout\n");
    //     bam_hdr_destroy(header);
    //     sam_close(in);
    //     return 1;
    // }

    // // Write header to output
    // if (sam_hdr_write(out, header) != 0) {
    //     fprintf(stderr, "Failed to write header\n");
    //     bam_hdr_destroy(header);
    //     sam_close(in);
    //     sam_close(out);
    //     return 1;
    // }

    // Initialize alignment record
    // bam1_t *aln = bam_init1();
    // if (!aln) {
    //     fprintf(stderr, "Failed to initialize alignment\n");
    //     bam_hdr_destroy(header);
    //     sam_close(in);
    //     // sam_close(out);
    //     return 1;
    // }
    // // 
    // FILE *fp = fopen("output.txt", "wb+");
    // if (fp == NULL) {
    //     perror("Failed to open file");
    //     return 1;
    // }
    // Writer* writer = create_writer(fp, header); // Assuming 1 is the file descriptor for stdout

    // if (writer == NULL)
    // {
    //     fprintf(stderr, "Failed to create writer\n");
    //     return 1;
    // }
    // // Read and write records
    // while (sam_read1(in, header, aln) >= 0) {
    //     if (write_bam_record(writer, aln) < 0) {
    //         fprintf(stderr, "Error writing alignment\n");
    //         break;
    //     }
    // }

    // close_writer(writer);


    // // Cleanup
    // bam_destroy1(aln);
    // bam_hdr_destroy(header);
    // sam_close(in);

    FILE *fp = fopen("output.txt", "rb");
    if (fp == NULL) {
        perror("Failed to open file");
        return 1;
    }
   
    bam1_t *aln = bam_init1();
    Reader* reader = make_reader(fp);
    for (int i = 0; i < reader->rec_num; i++) {

        // Convert to SAM format and print
        kstring_t str = {0, 0, NULL};
        read_record(reader, i, aln);
        sam_format1(reader->header, aln, &str);
        printf("%s\n", str.s);
        free(str.s);
        str.s = NULL;
        str.l = str.m = 0;
    }

    close_reader(reader);
    bam_destroy1(aln);
    
    return 0;
}
