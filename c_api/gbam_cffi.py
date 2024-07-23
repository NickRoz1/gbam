from cffi import FFI
import os

ffi = FFI()

# Define the C functions and types
ffi.cdef("""
void* get_gbam_reader(const char* gbam_path);
bam1_t get_bam_record(void* reader, uint64_u rec_num);
sam_hdr_t* get_header(void* reader);

void free_bam_record(bam1_t rec);
void free_gbam_reader(void* reader);
void free_header(sam_hdr_t* hdr);
""")

# Load the shared library
if os.name == "nt":
    lib = ffi.dlopen("../target/release/libgbam_tools_cffi.so")
else:
    raise OSError("Unsupported platform")


def get_reader(path):
    return lib.get_gbam_reader(path.encode('utf-8'))

def get_bam_record(reader, rec_num):
    return lib.get_bam_record(reader, rec_num)

def get_header(reader):
    return lib.get_header(reader)

# Example usage
if __name__ == "__main__":
    file_path = "/home/nickr/projects/gbam/test_data/one_record.gbam"
    reader = get_reader(file_path)
    
    if result == 0:
        print("File does not exist.")
        exit(1)

    hdr = get_header(reader)