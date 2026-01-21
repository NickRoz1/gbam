import argparse
import json
import struct
import mmap
import bisect
import lz4.block
import pysam
import zstandard as zstd
import brotli
from collections import OrderedDict
zstd_dctx = zstd.ZstdDecompressor()

# Hardcoded sizes for fixed fields (fallback if item_size is not in metadata)
FIXED_FIELD_SIZES = {
    "RefID": 4,
    "Pos": 4,
    "LName": 1,
    "Mapq": 1,
    "Bin": 2,
    "NCigar": 2,
    "Flags": 2,
    "SequenceLength": 4,
    "NextRefID": 4,
    "NextPos": 4,
    "TemplateLength": 4,
    "RawSeqLen": 4,
    "RawTagsLen": 4,
}

parser = argparse.ArgumentParser(description="Convert BAM to GBAM with column-wise compression.")
parser.add_argument("-i", "--input", required=True, help="Input BAM file path")
parser.add_argument("-o", "--output", required=True, help="Output GBAM file path")
args = parser.parse_args()

gbam_path = args.input
bam_input_path = args.output

# Open and memory-map the input GBAM file.
# gbam_path = "/Users/hasitha/Documents/biology/gbam/combined_compressed.gbam"
# gbam_path = "/Users/hasitha/Documents/biology/gbam/gbam_python/test_data/zstd_compressed.gbam"
f = open(gbam_path, "rb")
mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)

# Read header JSON from fixed-size region.
raw_header = mm[:1000].rstrip(b'\x00')
header_json = json.loads(raw_header.decode('utf-8'))
print("header_json:", json.dumps(header_json, indent=4))

# Read meta JSON (assumed to start at header_json["seekpos"])
mm.seek(header_json["seekpos"])
meta_json = json.loads(mm.read().decode('utf-8'))
print("meta_json:", json.dumps(meta_json, indent=4))

# --- Convert meta_json["field_to_meta"] into a dictionary keyed by field names ---
meta_field = meta_json["field_to_meta"]
field_to_meta = {}
if isinstance(meta_field, dict):
    field_to_meta = meta_field
elif isinstance(meta_field, list):
    for item in meta_field:
        if isinstance(item, dict):
            if "field" in item:
                field_to_meta[item["field"]] = item
            else:
                for key, value in item.items():
                    field_to_meta[key] = value
        elif isinstance(item, list) and len(item) == 2:
            field_to_meta[item[0]] = item[1]
        elif isinstance(item, str):
            print(f"Warning: Found a string in meta_json['field_to_meta']: {item}. Skipping.")
        else:
            print(f"Unexpected item type in meta_json['field_to_meta']: {type(item)}")
else:
    raise TypeError("meta_json['field_to_meta'] is neither a dict nor a list.")

if not field_to_meta:
    raise ValueError("No field metadata could be constructed from meta_json['field_to_meta'].")

# --- Helper: get blocks for a field ---
def get_blocks_for_field(field):
    meta = field_to_meta.get(field)
    if meta is None:
        raise KeyError(f"Field '{field}' not found in metadata (field_to_meta).")
    if isinstance(meta, dict):
        blocks = meta.get("blocks")
        if blocks is None:
            raise KeyError(f"No 'blocks' key in metadata for field '{field}'.")
        return blocks
    elif isinstance(meta, list):
        return meta
    else:
        raise TypeError(f"Unexpected type for metadata of field '{field}': {type(meta)}")

def encode_bam_tags(tag_list):
    import struct
    encoded = bytearray()
    for tag, value in tag_list:
        encoded.extend(tag.encode("utf-8"))
        if isinstance(value, int):
            encoded.extend(b'i')
            encoded.extend(struct.pack("<i", value))
        elif isinstance(value, float):
            encoded.extend(b'f')
            encoded.extend(struct.pack("<f", value))
        elif isinstance(value, str):
            encoded.extend(b'Z')
            encoded.extend(value.encode("utf-8") + b'\x00')
        else:
            raise TypeError(f"Unsupported tag value type: {tag} -> {type(value)}")
    return bytes(encoded)


# --- Optimized Column Classes ---

class OptimizedFixedColumn:
    def __init__(self, field):
        self.field = field
        self.blocks = get_blocks_for_field(field)
        if not self.blocks:
            raise ValueError(f"No block metadata for field: {field}")
        self.block_starts = []
        cum = 0
        for block in self.blocks:
            self.block_starts.append(cum)
            cum += block["numitems"]
        self.total_items = cum
        self.item_size = self.blocks[0].get("item_size", FIXED_FIELD_SIZES.get(field))
        self.decompressed = {}

    def get_decompressed_block(self, block_idx):
        if block_idx not in self.decompressed:
            block = self.blocks[block_idx]
            mm.seek(block["seekpos"])
            comp_buf = mm.read(block["block_size"])
            if block["block_size"] == block["uncompressed_size"]:
                data = comp_buf
            else:
                # data = lz4.block.decompress(comp_buf, uncompressed_size=block["uncompressed_size"])
                # data = zstd_dctx.decompress(comp_buf)
                # data = brotli.decompress(comp_buf)
                codec = block.get("codec", "lz4")
                if codec == "lz4":
                    data = lz4.block.decompress(comp_buf, uncompressed_size=block["uncompressed_size"])
                elif codec == "zstd":
                    zstd_dctx = zstd.ZstdDecompressor()
                    data = zstd_dctx.decompress(comp_buf)
                elif codec == "brotli":
                    data = brotli.decompress(comp_buf)
                else:
                    raise ValueError(f"Unsupported codec '{codec}' in block metadata for field '{self.field}'")
            self.decompressed[block_idx] = data
        return self.decompressed[block_idx]

    def get_item(self, i):
        block_idx = bisect.bisect_right(self.block_starts, i) - 1
        block_data = self.get_decompressed_block(block_idx)
        offset_in_block = (i - self.block_starts[block_idx]) * self.item_size
        return block_data[offset_in_block: offset_in_block + self.item_size]

class OptimizedVarColumn:
    def __init__(self, index_column, field):
        self.field = field
        self.blocks = get_blocks_for_field(field)
        if not self.blocks:
            raise ValueError(f"No block metadata for field: {field}")
        self.block_starts = []
        cum = 0
        for block in self.blocks:
            self.block_starts.append(cum)
            cum += block["numitems"]
        self.total_items = cum
        self.decompressed = {}
        self.index_column = index_column  # Provides length information.
        self.cum_cache = {}

    def get_decompressed_block(self, block_idx):
        if block_idx not in self.decompressed:
            block = self.blocks[block_idx]
            mm.seek(block["seekpos"])
            comp_buf = mm.read(block["block_size"])
            if block["block_size"] == block["uncompressed_size"]:
                data = comp_buf
            else:
                # data = lz4.block.decompress(comp_buf, uncompressed_size=block["uncompressed_size"])
                # data = zstd_dctx.decompress(comp_buf)
                # data = brotli.decompress(comp_buf)
                codec = block.get("codec", "lz4")
                if codec == "lz4":
                    data = lz4.block.decompress(comp_buf, uncompressed_size=block["uncompressed_size"])
                elif codec == "zstd":
                    zstd_dctx = zstd.ZstdDecompressor()
                    data = zstd_dctx.decompress(comp_buf)
                elif codec == "brotli":
                    data = brotli.decompress(comp_buf)
                else:
                    raise ValueError(f"Unsupported codec '{codec}' in block metadata for field '{self.field}'")
            self.decompressed[block_idx] = data
        return self.decompressed[block_idx]

    def get_cumulative(self, block_idx):
        if block_idx in self.cum_cache:
            return self.cum_cache[block_idx]
        start = self.block_starts[block_idx]
        num_items = self.blocks[block_idx]["numitems"]
        cum_list = [0] * num_items
        s = 0
        for j in range(start, start + num_items):
            val = self.index_column.get_item(j)[0]  # Each index item is one byte.
            s += val
            cum_list[j - start] = s
        self.cum_cache[block_idx] = cum_list
        return cum_list

    def get_item(self, i):
        block_idx = bisect.bisect_right(self.block_starts, i) - 1
        block_data = self.get_decompressed_block(block_idx)
        cum_list = self.get_cumulative(block_idx)
        rec_in_block = i - self.block_starts[block_idx]
        start_offset = 0 if rec_in_block == 0 else cum_list[rec_in_block - 1]
        end_offset = cum_list[rec_in_block]
        return block_data[start_offset:end_offset]

# --- Field Mapping ---
# Fixed fields are given as None; variable fields use the provided index.
fields_mapping = OrderedDict([
    ("RefID", None),
    ("Pos", None),
    ("LName", None),
    ("Mapq", None),
    ("Bin", None),
    ("NCigar", None),
    ("Flags", None),
    ("SequenceLength", None),
    ("NextRefID", None),
    ("NextPos", None),
    ("TemplateLength", None),
    ("RawSeqLen", None),
    ("RawTagsLen", None),
    ("ReadName", "LName"),    # Use LName as the length indicator.
    ("RawCigar", "NCigar"),   # Use NCigar as the length indicator.
    ("RawSequence", "RawSeqLen"),
    ("RawQual", "SequenceLength"),  # Use SequenceLength as the length indicator.
    ("RawTags", "RawTagsLen")
])

# Create column objects.
columns = {}
for field, index_field in fields_mapping.items():
    if index_field is None:
        columns[field] = OptimizedFixedColumn(field)
    else:
        index_col = OptimizedFixedColumn(index_field)
        columns[field] = OptimizedVarColumn(index_col, field)

# Determine the total number of records (using the "Flags" field metadata).
try:
    flags_blocks = get_blocks_for_field("Flags")
    records_num = sum(block["numitems"] for block in flags_blocks)
except KeyError as e:
    raise KeyError("The field 'Flags' is not present in metadata; cannot determine number of records.") from e
print(f"Total number of records: {records_num}")

# --- Rebuild the BAM header ---
header_text = bytes(meta_json["sam_header"])
L_text = len(header_text)
ref_seqs = meta_json["name_to_ref_id"]
n_ref = len(ref_seqs)

header_bin = bytearray()
header_bin.extend(b"BAM\x01")
header_bin.extend(struct.pack("<I", L_text))
header_bin.extend(header_text)
header_bin.extend(struct.pack("<I", n_ref))

for ref in ref_seqs:
    ref_name = ref[0].encode("utf-8") + b'\0'
    header_bin.extend(struct.pack("<I", len(ref_name)))
    header_bin.extend(ref_name)
    header_bin.extend(struct.pack("<I", ref[1]))

# --- Write the BAM file ---
# output_file_path = "/Users/hasitha/Documents/biology/gbam/reconverted_combined_compressed.bam"
with pysam.BGZFile(bam_input_path, "wb") as bam_file:
    bam_file.write(bytes(header_bin))
    for rec_i in range(records_num):
        arr = []
        # Fixed fields. LName NCigar RawSeqLen RawTagsLen
        arr.append(columns["RefID"].get_item(rec_i))
        arr.append(columns["Pos"].get_item(rec_i))
        read_name = columns["ReadName"].get_item(rec_i)
        arr.append(len(read_name).to_bytes(1, "little"))
        arr.append(columns["Mapq"].get_item(rec_i))
        arr.append(columns["Bin"].get_item(rec_i))
        raw_cigar = columns["RawCigar"].get_item(rec_i)
        arr.append(struct.pack('<H', len(raw_cigar) // 4))
        arr.append(columns["Flags"].get_item(rec_i))
        # Use the fixed "SequenceLength" field as l_seq.
        seq_len_bytes = columns["SequenceLength"].get_item(rec_i)
        seq_len = struct.unpack("<I", seq_len_bytes)[0]
        arr.append(struct.pack("<I", seq_len))
        arr.append(columns["NextRefID"].get_item(rec_i))
        arr.append(columns["NextPos"].get_item(rec_i))
        arr.append(columns["TemplateLength"].get_item(rec_i))
        # Variable fields.
        arr.append(read_name)
        arr.append(raw_cigar)
        raw_sequence = columns["RawSequence"].get_item(rec_i)
        arr.append(raw_sequence)
        raw_qual = columns["RawQual"].get_item(rec_i)
        # (Optional: check that raw_qual length matches seq_len)
        if len(raw_qual) != seq_len:
            print(f"Warning: record {rec_i} quality length {len(raw_qual)} != seq_len {seq_len}")
        arr.append(raw_qual)
        
        raw_tags_json = columns["RawTags"].get_item(rec_i)
        try:
            tag_list = json.loads(raw_tags_json.decode("utf-8"))
            tag_bytes = encode_bam_tags(tag_list)
            arr.append(tag_bytes)
        except Exception as e:
            print(f"Warning: Failed to decode tags for record {rec_i}: {e}")
            arr.append(b'')  # Add empty tag section on error

        total_len = sum(len(piece) for piece in arr)
        bam_file.write(struct.pack('<I', total_len))
        for piece in arr:
            bam_file.write(piece)
        if rec_i % 100000 == 0:
            print(f"Processed record {rec_i}")

print(f"BAM file successfully created: {bam_input_path}")
