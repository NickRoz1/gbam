use gbam_tools::reader::{parse_tmplt::ParsingTemplate, reader::Reader, record::GbamRecord};
use rust_htslib::htslib::{bam1_t, sam_hdr_destroy, sam_hdr_parse, sam_hdr_t, kstring_t};
use byteorder::{LittleEndian, ReadBytesExt};
use std::ffi::CStr;
use std::fs::File;
use std::path::Path;

#[no_mangle]
pub extern "C" fn get_gbam_reader(gbam_path: *const libc::c_char) -> *mut Reader {
    let gbam_path_str: String = unsafe {
        assert!(!gbam_path.is_null());
        CStr::from_ptr(gbam_path).to_string_lossy().into_owned()
    };

    let gbam_path = Path::new(&gbam_path_str);

    if !gbam_path.exists() {
        eprintln!(
            "The file doesnt exists or couldn't be found by Rust code: {}",
            gbam_path.to_str().unwrap()
        );
        return std::ptr::null_mut();
    }

    let file = File::open(gbam_path).unwrap();
    let mut tmplt = ParsingTemplate::new();
    tmplt.set_all();
    let reader = Reader::new(file, tmplt).unwrap();

    return Box::into_raw(Box::new(reader));
}

#[no_mangle]
pub extern "C" fn get_records_num(reader: *mut Reader) -> u64 {
    unsafe {
        (*reader).amount as u64
    }
}

#[no_mangle]
pub extern "C" fn get_bam_record(reader: *mut Reader, rec_num: u64) -> bam1_t {
    let mut rec = GbamRecord::default();
    unsafe {
        (*reader).fill_record(rec_num as usize, &mut rec);
    }

    return rec.get_hts_repr();
}

#[no_mangle]
pub extern "C" fn free_bam_record(rec: bam1_t) {
    unsafe {
        let _: Vec<u8> = Vec::from_raw_parts(rec.data, rec.l_data as usize, rec.m_data as usize);
    }
}

#[no_mangle]
pub extern "C" fn free_gbam_reader(reader: *mut Reader) {
    unsafe {
        drop(Box::<Reader>::from_raw(reader));
    }
}

#[no_mangle]
pub extern "C" fn get_header(reader: *mut Reader) -> *mut sam_hdr_t {
    unsafe {
        let header_len = (&(*reader).file_meta.get_sam_header()[..std::mem::size_of::<u32>()]).read_u32::<LittleEndian>().unwrap() as usize;
        let header_bytes = (*reader).file_meta.get_sam_header()[std::mem::size_of::<u32>()..std::mem::size_of::<u32>()+header_len].to_owned();

        sam_hdr_parse(header_len as u64, header_bytes.as_ptr() as *const i8)
    }
}


#[no_mangle]
pub extern "C" fn free_header(hdr: *mut sam_hdr_t) {
    unsafe{
        sam_hdr_destroy(hdr);
    }   
}

// int sam_parse1(kstring_t *s, sam_hdr_t *h, bam1_t *b) HTS_RESULT_USED;
// Make a python test
// Need to handle the header
// Verify in Python that the C API works fine and there are no leaks
// Integrate this into SAMTOOLS view code and verify it works
// Integrate this into SAMTOOLS depth and SAMTOOLS markdup


// The serializer from bincode does bullshit, it's just for saving not for exchangin data with samtools
// Pass the data to samtools directly
// The problem to figure out is the header, it has to be encoded and passed to samtools somehow too


// Create a header from existing text.
//
// @param l_text    Length of text
// @param text      Header text
// @return A populated sam_hdr_t structure on success; NULL on failure.
// @note The text field of the returned header will be NULL, and the l_text
// field will be zero.
//
// The sam_hdr_t struct returned by a successful call should be freed
// via sam_hdr_destroy() when it is no longer needed.
//
// HTSLIB_EXPORT
// sam_hdr_t *sam_hdr_parse(size_t l_text, const char *text);
// use rust_htslib::htslib::{bam_write1, sam_hdr_t, sam_hdr_parse};
