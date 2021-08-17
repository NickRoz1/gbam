use std::sync::Arc;
use std::{borrow::Borrow, fs::File, rc::Rc};

use bam_tools::record::fields::{
    field_type, var_size_field_to_index, FieldType, Fields, FIELDS_NUM,
};
use memmap2::Mmap;

use crate::meta::{FileInfo, FileMeta, FILE_INFO_SIZE};
use crate::writer::calc_crc_for_meta_bytes;

use super::{
    column::{Column, FixedColumn, Inner, VariableColumn},
    parse_tmplt::ParsingTemplate,
    record::GbamRecord,
    records::Records,
};

pub struct Reader {
    // Instead of hashmap. Empty columns will contain None.
    pub columns: Vec<Option<Box<dyn Column + Send>>>,
    pub parsing_template: ParsingTemplate,
    pub rec_num: usize,
    pub file_meta: Arc<FileMeta>,
    // Kept so File won't drop while used by mmap.
    inner: Box<File>,
}

impl Reader {
    pub fn new(inner: File, parsing_template: ParsingTemplate) -> std::io::Result<Self> {
        let inner = Box::new(inner);
        let mmap = Arc::new(unsafe { Mmap::map(inner.borrow())? });
        let file_meta = Arc::new(verify_and_parse_meta(&mmap)?);
        let rec_num = file_meta
            .view_blocks(&Fields::RefID)
            .iter()
            .fold(0, |acc, x| acc + x.numitems) as usize;
        Ok(Self {
            columns: init_columns(&mmap, &parsing_template, &file_meta),
            parsing_template,
            file_meta,
            rec_num,
            inner,
        })
    }

    pub fn fill_record(&mut self, rec_num: usize, rec: &mut GbamRecord) {
        for &field in self.parsing_template.get_active_data_fields_iter() {
            self.columns[field as usize]
                .as_mut()
                .unwrap()
                .fill_record_field(rec_num, rec);
        }
    }

    /// Get iterator over all GBAM records (according to parsing template).
    pub fn records(&mut self) -> Records {
        Records::new(self)
    }
}

fn init_columns(
    mmap: &Arc<Mmap>,
    parse_template: &ParsingTemplate,
    meta: &Arc<FileMeta>,
) -> Vec<Option<Box<dyn Column + Send>>> {
    let mut res = Vec::new();
    (0..FIELDS_NUM).for_each(|_| res.push(None));
    for &field in parse_template.get_active_fields_iter() {
        res[field as usize] = Some(init_col(field, mmap, meta));
    }
    res
}

fn init_col(field: Fields, mmap: &Arc<Mmap>, meta: &Arc<FileMeta>) -> Box<dyn Column + Send> {
    let inner = Inner::new(meta.clone(), field, mmap.clone());
    match field_type(&field) {
        FieldType::FixedSized => Box::new(FixedColumn::new(inner)),
        FieldType::VariableSized => {
            let idx_field = var_size_field_to_index(&field);
            let idx_inner = Inner::new(meta.clone(), idx_field, mmap.clone());
            let idx_col = FixedColumn::new(idx_inner);
            Box::new(VariableColumn::new(inner, idx_col))
        }
    }
}

fn verify_and_parse_meta(mmap: &Mmap) -> std::io::Result<FileMeta> {
    let file_info_bytes = &mmap[0..FILE_INFO_SIZE];
    let file_info = FileInfo::from(&file_info_bytes[..]);
    // Read file meta
    let buf = &mmap[file_info.seekpos as usize..];
    if calc_crc_for_meta_bytes(&buf[..]) != file_info.crc32 {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "Metadata JSON was damaged.",
        ));
    }
    let file_meta_json_str = String::from_utf8(buf.to_owned()).unwrap();
    Ok(serde_json::from_str(&file_meta_json_str).expect("File meta json string was damaged."))
}
