#include <htslib/sam.h>
#include <stdlib.h>
#include <unistd.h>
#include "defs.h"

int write_gbam(char* file_path){
    // Open the file for writing
    FILE *fp = fopen(file_path, "wb");
    if (!fp) {
        perror("Failed to open file for writing");
        return;
    }

    samFile *in = sam_open("-", "r");
    if (!in) {
        fprintf(stderr, "Failed to open stdin\n");
        return 1;
    }

    if (hts_set_threads(in, 8) != 0) {
        sam_close(fp);
        return NULL;
    }

    bam_hdr_t *header = sam_hdr_read(in);
    if (!header) {
        fprintf(stderr, "Failed to read header\n");
        sam_close(in);
        return 1;
    }

    bam1_t *aln = bam_init1();
    if (!aln) {
        fprintf(stderr, "Failed to initialize alignment\n");
        bam_hdr_destroy(header);
        sam_close(in);
        // sam_close(out);
        return 1;
    }

    Writer* writer = create_writer(fp, header); // Assuming 1 is the file descriptor for stdout
    if(writer == NULL) {
        fprintf(stderr, "Failed to create writer\n");
        bam_destroy1(aln);
        bam_hdr_destroy(header);
        sam_close(in);
        fclose(fp);
        return 1;
    }

    while (sam_read1(in, header, aln) >= 0) {
        if (write_bam_record(writer, aln) < 0) {
            fprintf(stderr, "Error writing alignment\n");
            bam_destroy1(aln);
            bam_hdr_destroy(header);
            sam_close(in);
            fclose(fp);
            return 1;
        }
    }

    close_writer(writer);

    // Cleanup
    bam_destroy1(aln);
    bam_hdr_destroy(header);
    sam_close(in);
    fclose(fp);

    return 0;
}

int read_gbam(char* file_path){
    // Open the file for reading
    FILE *fp = fopen(file_path, "rb");
    if (!fp) {
        perror("Failed to open file for reading");
        return 1;
    }

    Reader* reader = make_reader(fp);
    if (reader == NULL) {
        fprintf(stderr, "Failed to create reader\n");
        fclose(fp);
        return 1;
    }

    bam1_t *aln = bam_init1();
    if (!aln) {
        fprintf(stderr, "Failed to initialize alignment\n");
        close_reader(reader);
        fclose(fp);
        return 1;
    }

    for (int i = 0; i < reader->rec_num; i++) {
        read_record(reader, i, aln);

        // Convert to SAM format and print
        kstring_t str = {0, 0, NULL};
        sam_format1(reader->header, aln, &str);
        printf("%s\n", str.s);
        free(str.s);
        str.s = NULL;
        str.l = str.m = 0;
    }

    close_reader(reader);
    bam_destroy1(aln);
    fclose(fp);

    return 0;
}

int main(int argc, char *argv[]) {

    if (argc < 2) {
        fprintf(stderr, "Usage: %s <file-path>.gbam\nYou can also pipe sam/bam into the tool and then it will instead create a gbam file.\n", argv[0]);
        return 1;
    }

    char *file_path = argv[1];

    if (!isatty(fileno(stdin))) {
        return write_gbam(file_path);
    }

    return read_gbam(file_path);
}
