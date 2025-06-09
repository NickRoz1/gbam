#include <htslib/sam.h>
#include <stdlib.h>
#include <unistd.h>
#include "defs.h"
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdio.h>
#include <sys/stat.h>

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

    int fd = open(file_path, O_RDONLY);
    if (fd == -1) {
        perror("Error opening file");
        return 1;
    }

    struct stat sb;
    fstat(fd, &sb);
    size_t file_size = sb.st_size;
    
    // Map file into memory (read-only)
    char *mapped = mmap(NULL, file_size, PROT_READ, MAP_SHARED, fd, 0);
    if (mapped == MAP_FAILED) {
        perror("Error mapping file");
        close(fd);
        return 1;
    }

    Reader* reader = make_reader(mapped);

    bam1_t *aln = bam_init1();

    char stdout_buffer[65536];
    flockfile(stdout);

    setvbuf(stdout, stdout_buffer, _IOFBF, sizeof(stdout_buffer));
    kstring_t str = {0, 0, NULL};

    for (int i = 0; i < reader->rec_num; i++) {
        read_record(reader, i, aln);

        // Convert to SAM format and print
        sam_format1(reader->header, aln, &str);
        printf("%s\n", str.s);
        str.l = 0;
    }
    free(str.s);


    funlockfile(stdout);
    close_reader(reader);
    bam_destroy1(aln);
    munmap(mapped, file_size);
    close(fd);

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
