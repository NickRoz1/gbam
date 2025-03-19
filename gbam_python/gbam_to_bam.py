import json
from collections import OrderedDict
import lz4.block as lb
import lz4.frame
import pysam
import struct
import sys
import functools 

f = open("/Users/hasitha/Documents/biology/gbam/gbam_python/output.gbam", "rb")
header = (f.read(1000).rstrip(b'\x00'))
header_json = json.loads(header.decode('utf-8'))
print("header_json:", json.dumps(header_json, indent=4))

f.seek(header_json["seekpos"])
meta_json = json.loads(f.read().decode('utf-8'))
print("meta_json:", json.dumps(meta_json, indent=4))

fields_mapping = OrderedDict()
fields_mapping["RefID"] = None
fields_mapping["Pos"] = None
fields_mapping["LName"] = None
fields_mapping["Mapq"] = None
fields_mapping["Bin"] = None
fields_mapping["NCigar"] = None
fields_mapping["Flags"] = None
fields_mapping["SequenceLength"] = None
fields_mapping["NextRefID"] = None
fields_mapping["NextPos"] = None
fields_mapping["TemplateLength"] = None
fields_mapping["RawSeqLen"] = None
fields_mapping["RawTagsLen"] = None
fields_mapping["ReadName"] = "LName"
fields_mapping["RawCigar"] = "NCigar"
fields_mapping["RawSequence"] = "RawSeqLen"
fields_mapping["RawQual"] = "SequenceLength"
fields_mapping["RawTags"] = "RawTagsLen"

@functools.cache
def get_decompressed(seekpos, block_size, uncompressed_size):
    print(f"Decompressing block at seekpos={seekpos}, block_size={block_size}, uncompressed_size={uncompressed_size}")
    
    f.seek(seekpos)
    compressed_buf = f.read(block_size)
    
    print(f"File position after reading: {f.tell()}, buffer length: {len(compressed_buf)}")
    print(f"Compressed buffer (first 20 bytes): {compressed_buf[:20].hex()}")

    if block_size == uncompressed_size:
        print("Skipping decompression, using raw data.")
        return compressed_buf  # Directly return raw data if it's not compressed
    
    try:
        decompressed_data = lb.decompress(compressed_buf, uncompressed_size)
        print(f"Decompressed block size: {len(decompressed_data)}")
        return decompressed_data
    except Exception as e:
        print(f"ERROR: LZ4 decompression failed at seekpos={seekpos}")
        print(f"Exception: {e}")
        sys.exit(1)

def get_block(i, meta):
    blocks = meta["blocks"]
    loaded = blocks[0]["numitems"]
    cur = 0
    while loaded <= i:
        cur += 1
        loaded += blocks[cur]["numitems"]
    
    # print(f"Fetching block {cur}, numitems={blocks[cur]['numitems']}")
    dec = get_decompressed(blocks[cur]["seekpos"], blocks[cur]["block_size"], blocks[cur]["uncompressed_size"])
    assert len(dec) == blocks[cur]["uncompressed_size"], "Decompressed size mismatch!"
    return (loaded - blocks[cur]["numitems"], dec)

# Define fixed field sizes (in bytes) for the fixed fields.
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

class FixedColumn():
    def __init__(self, field):
        self.field = field
        # Wrap the list in a dictionary under "blocks"
        self.meta = {"blocks": meta_json["field_to_meta"][field]}

    def get_item(self, i):
        loaded, uncompressed_buf = get_block(i, self.meta)
        i -= loaded
        # Try to get item_size from metadata.
        item_size = self.meta["blocks"][0].get("item_size")
        if item_size is None:
            # Fallback: use the hardcoded fixed field size.
            item_size = FIXED_FIELD_SIZES[self.field]
        offset = i * item_size
        return uncompressed_buf[offset:(offset + item_size)]

class VarColumn():
    def __init__(self, fixed_index_column, field):
        self.meta = {"blocks": meta_json["field_to_meta"][field]}
        self.index = fixed_index_column

    def get_item(self, i):
        loaded, uncompressed_buf = get_block(i, self.meta)
        # Compute cumulative offset by summing one-byte values from index
        offset = 0
        for j in range(loaded, i + 1):
            offset += self.index.get_item(j)[0]  # each returns a 1-byte value
        # Also compute the previous cumulative offset:
        prev_offset = 0
        for j in range(loaded, i):
            prev_offset += self.index.get_item(j)[0]
        return uncompressed_buf[prev_offset:offset]

    def get_block(i, meta):
        blocks = meta["blocks"]
        loaded = blocks[0]["numitems"]
        cur = 0
        while loaded <= i:
            cur += 1
            loaded += blocks[cur]["numitems"]
        
        print(f"Fetching block {cur}, numitems={blocks[cur]['numitems']}")
        dec = get_decompressed(blocks[cur]["seekpos"], blocks[cur]["block_size"], blocks[cur]["uncompressed_size"])
        assert len(dec) == blocks[cur]["uncompressed_size"], "Decompressed size mismatch!"
        return (loaded - blocks[cur]["numitems"], dec)

records_num = sum(block["numitems"] for block in meta_json["field_to_meta"]["Flags"])
print(f"Total number of records: {records_num}")

columns = {}
for key, value in fields_mapping.items():
    if value is None:
        columns[key] = FixedColumn(key)
    else:
        index = FixedColumn(value)
        columns[key] = VarColumn(index, key)

# Convert the SAM header from meta_json (assuming it's stored as a list of ints)
header_text = bytes(meta_json["sam_header"])
L_text = len(header_text)

# Get reference sequence info from the metadata.
# meta_json["name_to_ref_id"] is assumed to be a list of [ref_name, ref_length] pairs.
ref_seqs = meta_json["name_to_ref_id"]
n_ref = len(ref_seqs)

# Build the binary header.
header_bin = bytearray()
header_bin.extend(b"BAM\x01")
header_bin.extend(struct.pack("<I", L_text))  # L_text: length of SAM header text
header_bin.extend(header_text)
header_bin.extend(struct.pack("<I", n_ref))  # Number of reference sequences

for ref in ref_seqs:
    # ref[0] is the reference name; add a trailing null byte.
    ref_name = ref[0].encode("utf-8") + b'\0'
    # The header requires the length of the reference name (including the null)
    header_bin.extend(struct.pack("<I", len(ref_name)))
    header_bin.extend(ref_name)
    # Write the reference length as a 4-byte little-endian integer.
    header_bin.extend(struct.pack("<I", ref[1]))

# Now write the header to the output BAM file.
output_file_path = "/Users/hasitha/Documents/biology/gbam/output.bam"
# with open(output_file_path, "wb") as bam_file:
#     bam_file.write(header_bin)
    
#     # Write the rest of the records as before.
#     for rec_i in range(0, records_num):
#         arr = []
#         arr.append(columns["RefID"].get_item(rec_i))
#         arr.append(columns["Pos"].get_item(rec_i))
#         arr.append(len(columns["ReadName"].get_item(rec_i)).to_bytes(1, "little"))
#         arr.append(columns["Mapq"].get_item(rec_i))
#         arr.append(columns["Bin"].get_item(rec_i))
#         arr.append(struct.pack('<H', len(columns["RawCigar"].get_item(rec_i)) // 4))
#         arr.append(columns["Flags"].get_item(rec_i))
#         arr.append(struct.pack('<L', len(columns["RawQual"].get_item(rec_i))))
#         arr.append(columns["NextRefID"].get_item(rec_i))
#         arr.append(columns["NextPos"].get_item(rec_i))
#         arr.append(columns["TemplateLength"].get_item(rec_i))
#         arr.append(columns["ReadName"].get_item(rec_i))
#         arr.append(columns["RawCigar"].get_item(rec_i))
#         arr.append(columns["RawSequence"].get_item(rec_i))
#         arr.append(columns["RawQual"].get_item(rec_i))
#         # arr.append(columns["RawTags"].get_item(rec_i))
        
#         total_len = sum(len(a) for a in arr)
#         arr.insert(0, struct.pack('<L', total_len))
        
#         for a in arr:
#             bam_file.write(a)
            
#         if rec_i == 100:
#             break
#     bgzf_eof = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000")
#     bam_file.write(bgzf_eof)
with pysam.BGZFile(output_file_path, "wb") as bam_file:
    bam_file.write(bytes(header_bin))
    for rec_i in range(0, records_num):
        arr = []
        arr.append(columns["RefID"].get_item(rec_i))
        arr.append(columns["Pos"].get_item(rec_i))
        arr.append(len(columns["ReadName"].get_item(rec_i)).to_bytes(1, "little"))
        arr.append(columns["Mapq"].get_item(rec_i))
        arr.append(columns["Bin"].get_item(rec_i))
        arr.append(struct.pack('<H', len(columns["RawCigar"].get_item(rec_i)) // 4))
        arr.append(columns["Flags"].get_item(rec_i))
        arr.append(struct.pack('<L', len(columns["RawQual"].get_item(rec_i))))
        arr.append(columns["NextRefID"].get_item(rec_i))
        arr.append(columns["NextPos"].get_item(rec_i))
        arr.append(columns["TemplateLength"].get_item(rec_i))
        arr.append(columns["ReadName"].get_item(rec_i))
        arr.append(columns["RawCigar"].get_item(rec_i))
        arr.append(columns["RawSequence"].get_item(rec_i))
        arr.append(columns["RawQual"].get_item(rec_i))
        # If you have aux data in RawTags, append it here.
        total_len = sum(len(a) for a in arr)
        arr.insert(0, struct.pack('<L', total_len))
        for a in arr:
            bam_file.write(a)
            
        if rec_i == 100:
            break
print(f"BAM file successfully created: {output_file_path}")
