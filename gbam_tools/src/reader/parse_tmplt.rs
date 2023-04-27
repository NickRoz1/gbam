#[cfg(feature = "python-ffi")]
use pyo3::prelude::*;

use bam_tools::record::fields::{
    field_type, is_data_field, var_size_field_to_index, FieldType, Fields, DATA_FIELDS_NUM,
    FIELDS_NUM,
};

/// This struct regulates what fields are getting parsed from GBAM file.
#[cfg_attr(feature = "python-ffi", pyclass)]
#[derive(Clone, Debug)]
pub struct ParsingTemplate {
    inner: Vec<Option<Fields>>,
    // Cache.
    active_data_fields: Vec<Fields>,
}

impl ParsingTemplate {
    /// Create new parsing templates with all fields set to false
    pub fn new() -> Self {
        Self {
            inner: ((0..FIELDS_NUM).map(|_| None).collect()),
            active_data_fields: Vec::new(),
        }
    }

    /// Create new parsing templates with passed fields set to active
    pub fn new_with(fields_to_set: &[Fields]) -> Self {
        let mut empty = Self::new();
        for field in fields_to_set {
            empty.set(field, true);
        }
        empty
    }
    /// Set field value
    pub fn set(&mut self, field: &Fields, val: bool) {
        match field_type(field) {
            FieldType::FixedSized => {
                self.inner[*field as usize] = Self::bool_to_val(field, val);
            }
            FieldType::VariableSized => {
                self.inner[*field as usize] = Self::bool_to_val(field, val);
                self.inner[var_size_field_to_index(field) as usize] =
                    Self::bool_to_val(&var_size_field_to_index(field), val);
            }
        };
        self.set_active();
    }

    fn bool_to_val(field: &Fields, val: bool) -> Option<Fields> {
        match val {
            true => Some(*field),
            false => None,
        }
    }
    /// Get iterator over fields currently requested for parsing
    #[allow(clippy::needless_lifetimes)]
    pub fn get_active_fields_iter<'a>(&'a self) -> impl Iterator<Item = &'a Fields> {
        self.inner
            .iter()
            .filter(|x| x.is_some())
            .map(|x| x.as_ref().unwrap())
    }

    /// Get iterator over data fields (no index fields) currently requested for parsing
    #[allow(clippy::needless_lifetimes)]
    pub fn get_active_data_fields_iter<'a>(&'a self) -> impl Iterator<Item = &'a Fields> {
        self.active_data_fields.iter()
    }

    /// This method exists because list of active fields is stored separately.
    fn set_active(&mut self) {
        self.active_data_fields = self
            .inner
            .iter()
            .filter(|x| x.is_some() && is_data_field(&x.unwrap()))
            .map(|x| x.unwrap())
            .collect();
    }

    /// Get fields currently requested for parsing
    pub fn get_active_fields(&self) -> Vec<Fields> {
        self.get_active_fields_iter()
            .map(|field| field.to_owned())
            .collect::<Vec<_>>()
    }
    /// Set all fields to active state
    pub fn set_all(&mut self) {
        for (field, val) in Fields::iterator().zip(self.inner.iter_mut()) {
            if val.is_none() {
                *val = Some(*field);
            }
        }
        self.set_active();
    }
    /// Set all fields to disabled state
    pub fn clear(&mut self) {
        self.inner
            .iter_mut()
            .filter(|x| x.is_some())
            .for_each(|e| *e = None);
        self.set_active();
    }

    pub fn check_if_active(&self, fields: &[Fields]) -> bool {
        for &field in fields.iter() {
            if self.inner[field as usize].is_none() {
                return false;
            }
        }
        true
    }
}

impl Default for ParsingTemplate {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(feature = "python-ffi")]
#[pymethods]
impl ParsingTemplate {
    #[new]
    /// Create new ParsingTemplate with predefined values
    pub fn new_from_python(fields: Vec<bool>) -> Self {
        assert_eq!(fields.len(), DATA_FIELDS_NUM);
        let mut tmplt = ParsingTemplate::new();
        for (field, val) in Fields::iterator().zip(fields) {
            tmplt.set(field, val);
        }
        tmplt
    }
}
