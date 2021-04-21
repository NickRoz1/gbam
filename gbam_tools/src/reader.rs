use super::{Fields, FIELDS_NUM, GBAM_MAGIC};
use std::collections::HashMap;
use std::io::Read;
use std::rc::Rc;

/// GBAM record only includes queried fields.
pub struct GbamRecord {
    // Used to determine which fields are currently available.
    offsets: HashMap<Fields, usize>,
    // The loaded fields.
    data: Vec<u8>,
}

impl GbamRecord {
    pub fn new(offsets: HashMap<Fields, usize>, data: Vec<u8>) -> Self {
        GbamRecord {
            offsets: offsets,
            data: data,
        }
    }

    /// Checks whether field is parsed.
    /// Panics on fail.
    fn is_parsed(&self, field: Fields) {
        assert!(self.offsets.contains_key(&field));
    }

    /// Read refID
    pub fn read_refID(&self) {
        self.is_parsed(Fields::RefID);
    }

    /// Read pos
    pub fn read_pos(&mut self) {
        self.is_parsed(Fields::Pos);
    }

    /// Read length of read name
    pub fn read_l_read_name(&mut self) {
        self.is_parsed(Fields::LName);
    }

    /// Read mapq
    pub fn read_mapq(&mut self) {
        self.is_parsed(Fields::Mapq);
    }

    /// Read bin
    pub fn read_bin(&mut self) {
        self.is_parsed(Fields::Bin);
    }

    /// Read number of cigar operations
    pub fn read_ncigar(&mut self) {
        self.is_parsed(Fields::NCigar);
    }

    /// Read flags
    pub fn read_flags(&mut self) {
        self.is_parsed(Fields::Flags);
    }

    /// Read length of sequence
    pub fn read_l_seq(&mut self) {
        self.is_parsed(Fields::SequenceLength);
    }

    /// Read next refID
    pub fn read_next_refID(&mut self) {
        self.is_parsed(Fields::NextRefID);
    }

    /// Read next position
    pub fn read_next_pos(&mut self) {
        self.is_parsed(Fields::NextPos);
    }

    /// Read length of template
    pub fn read_l_template(&mut self) {
        self.is_parsed(Fields::TemplateLength);
    }

    /// Read read name
    pub fn read_read_name(&mut self) {
        self.is_parsed(Fields::Bin);
        self.0[Fields::ReadName as usize] = true;
        self.0[Fields::LName as usize] = true;
    }

    /// Read raw cigar
    pub fn read_raw_cigar(&mut self) {
        self.is_parsed(Fields::Bin);
        self.0[Fields::RawCigar as usize] = true;
        self.0[Fields::NCigar as usize] = true;
    }

    /// Read raw sequence
    pub fn read_raw_sequence(&mut self) {
        self.is_parsed(Fields::Bin);
        self.0[Fields::RawSequence as usize] = true;
        self.0[Fields::SequenceLength as usize] = true;
    }

    /// Read base qualities
    pub fn read_raw_qual(&mut self) {
        self.is_parsed(Fields::Bin);
        self.0[Fields::RawQual as usize] = true;
        self.0[Fields::SequenceLength as usize] = true;
    }

    /// Read tags
    pub fn read_raw_tags(&mut self) {
        self.is_parsed(Fields::Bin);
        self.0[Fields::RawTags as usize] = true;
        self.0[Fields::RawTagsLen as usize] = true;
    }

    /// Read length of tags
    pub fn read_l_tags(&mut self) {
        self.is_parsed(Fields::Bin);
        self.0[Fields::RawTagsLen as usize] = true;
    }
}
/// GBAM reader is capable of parsing whole rowgroups (assembling the splitted
/// fields of BAM records back into full record), or quickly parsing a few
/// columns (BAM fields ignoring others).
pub struct Reader {
    reading_mode: ReadingMode,
    buffer: Vec<GbamRecord>,
}

impl Reader {
    pub fn new(reading_mode: ReadingMode) -> Self {
        Reader {
            reading_mode: reading_mode,
        }
    }
}

/// This structure is used to customize reading mode for GBAM reader. By calling
/// methods of this struct one can customize which fields will be parsed.
pub struct ReadingMode(Vec<bool>);

impl ReadingMode {
    pub fn new() -> Self {
        Self(vec![false; FIELDS_NUM])
    }
    /// Assembles all BAM records from rowgroup batch back.
    pub fn read_all(&mut self) {
        for setting in self.0.iter_mut() {
            *setting = true;
        }
    }

    /// Read refID
    pub fn read_refID(&mut self) {
        self.0[Fields::RefID as usize] = true;
    }

    /// Read pos
    pub fn read_pos(&mut self) {
        self.0[Fields::Pos as usize] = true;
    }

    /// Read length of read name
    pub fn read_l_read_name(&mut self) {
        self.0[Fields::LName as usize] = true;
    }

    /// Read mapq
    pub fn read_mapq(&mut self) {
        self.0[Fields::Mapq as usize] = true;
    }

    /// Read bin
    pub fn read_bin(&mut self) {
        self.0[Fields::Bin as usize] = true;
    }

    /// Read number of cigar operations
    pub fn read_ncigar(&mut self) {
        self.0[Fields::NCigar as usize] = true;
    }

    /// Read flags
    pub fn read_flags(&mut self) {
        self.0[Fields::Flags as usize] = true;
    }

    /// Read length of sequence
    pub fn read_l_seq(&mut self) {
        self.0[Fields::SequenceLength as usize] = true;
    }

    /// Read next refID
    pub fn read_next_refID(&mut self) {
        self.0[Fields::NextRefID as usize] = true;
    }

    /// Read next position
    pub fn read_next_pos(&mut self) {
        self.0[Fields::NextPos as usize] = true;
    }

    /// Read length of template
    pub fn read_l_template(&mut self) {
        self.0[Fields::TemplateLength as usize] = true;
    }

    /// Read read name
    pub fn read_read_name(&mut self) {
        self.0[Fields::ReadName as usize] = true;
        self.0[Fields::LName as usize] = true;
    }

    /// Read raw cigar
    pub fn read_raw_cigar(&mut self) {
        self.0[Fields::RawCigar as usize] = true;
        self.0[Fields::NCigar as usize] = true;
    }

    /// Read raw sequence
    pub fn read_raw_sequence(&mut self) {
        self.0[Fields::RawSequence as usize] = true;
        self.0[Fields::SequenceLength as usize] = true;
    }

    /// Read base qualities
    pub fn read_raw_qual(&mut self) {
        self.0[Fields::RawQual as usize] = true;
        self.0[Fields::SequenceLength as usize] = true;
    }

    /// Read tags
    pub fn read_raw_tags(&mut self) {
        self.0[Fields::RawTags as usize] = true;
        self.0[Fields::RawTagsLen as usize] = true;
    }

    /// Read length of tags
    pub fn read_l_tags(&mut self) {
        self.0[Fields::RawTagsLen as usize] = true;
    }
}
