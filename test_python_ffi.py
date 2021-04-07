from cffi import FFI
ffi = FFI()
ffi.cdef("""
    void bam_to_gbam_python(char*, char*);
""")

C = ffi.dlopen("target/debug/libgbam_tools.so")

b_string1 = "Hi".encode('utf-8')
b_string2 = "It's working".encode('utf-8')

C.bam_to_gbam_python(b_string1, b_string2)