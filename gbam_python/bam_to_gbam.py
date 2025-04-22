import argparse
import json
import struct
import lz4.block
import pysam
import zstandard as zstd
import brotli
from collections import OrderedDict

# Compression codec options
COMPRESSION_LZ4 = "lz4"
COMPRESSION_ZSTD = "zstd"
COMPRESSION_BROTLI = "brotli"

# Choose compression type per column
column_compression = {
    "RefID": COMPRESSION_LZ4,
    "Pos": COMPRESSION_LZ4,
    "LName": COMPRESSION_LZ4,
    "Mapq": COMPRESSION_LZ4,
    "Bin": COMPRESSION_LZ4,
    "NCigar": COMPRESSION_LZ4,
    "Flags": COMPRESSION_LZ4,
    "SequenceLength": COMPRESSION_LZ4,
    "NextRefID": COMPRESSION_LZ4,
    "NextPos": COMPRESSION_LZ4,
    "TemplateLength": COMPRESSION_LZ4,
    "RawSeqLen": COMPRESSION_ZSTD,
    "RawTagsLen": COMPRESSION_ZSTD,
    "ReadName": COMPRESSION_BROTLI,
    "RawCigar": COMPRESSION_BROTLI,
    "RawSequence": COMPRESSION_ZSTD,
    "RawQual": COMPRESSION_ZSTD,
    "RawTags": COMPRESSION_BROTLI
}

# --- Configuration ---
BLOCK_SIZE_LIMIT = 1024 * 1024  # example block size; adjust as needed

# Define the fields mapping
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

# Initialize a dict to hold a bytearray buffer for each field.
columns = {field: bytearray() for field in fields_mapping}

# For each field, we will also record a list of block metadata.
field_meta = {field: [] for field in fields_mapping}

parser = argparse.ArgumentParser(description="Convert BAM to GBAM with column-wise compression.")
parser.add_argument("-i", "--input", required=True, help="Input BAM file path")
parser.add_argument("-o", "--output", required=True, help="Output GBAM file path")
args = parser.parse_args()

bam_input_path = args.input
gbam_output_path = args.output

# --- Helper: flush a column buffer if it exceeds a limit ---
def flush_column(field, out_f):
    buf = columns[field]
    if not buf:
        return
    uncompressed_data = bytes(buf)
    # (If desired,  may choose to skip compression if block size equals uncompressed size.)
    if len(uncompressed_data) > 0:
        print("COMPRESSED")
    # compressed_data = lz4.block.compress(uncompressed_data, store_size=False) if len(uncompressed_data) > 0 else uncompressed_data
    # cctx = zstd.ZstdCompressor(level=3)  # level 1–22
    # compressed_data = cctx.compress(uncompressed_data) if len(uncompressed_data) > 0 else uncompressed_data
    # compressed_data = brotli.compress(uncompressed_data, quality=5)
    
    codec = column_compression.get(field, COMPRESSION_LZ4)

    if codec == COMPRESSION_LZ4:
        compressed_data = lz4.block.compress(uncompressed_data, store_size=False)
    elif codec == COMPRESSION_ZSTD:
        cctx = zstd.ZstdCompressor(level=3)
        compressed_data = cctx.compress(uncompressed_data)
    elif codec == COMPRESSION_BROTLI:
        compressed_data = brotli.compress(uncompressed_data, quality=5)
    else:
        raise ValueError(f"Unknown compression codec: {codec} for field: {field}")

    seekpos = out_f.tell()
    out_f.write(compressed_data)
    meta = {
        "seekpos": seekpos,
        "numitems": record_count,  # total records written so far (for this column)
        "block_size": len(compressed_data),
        "uncompressed_size": len(uncompressed_data),
        "codec": codec,
        "stats": None  #  can add statistics if needed
    }
    field_meta[field].append(meta)
    # Clear the column buffer for new data.
    columns[field].clear()

# --- Open the BAM file for reading ---
bam_file = pysam.AlignmentFile(bam_input_path, "rb")

# Extract SAM header (as bytes) and reference sequences.
sam_header_bytes = bam_file.text.encode('utf-8')
ref_seqs = [(sq['SN'], sq['LN']) for sq in bam_file.header.get('SQ', [])]

record_count = 0

# --- Process each record ---
for rec in bam_file.fetch(until_eof=True):
    record_count += 1

    # Fixed fields:
    # Note: We pack integers in little-endian format.
    refid = struct.pack('<i', rec.reference_id)
    pos = struct.pack('<i', rec.reference_start)
    # LName: 1-byte length of the read name
    read_name = rec.query_name.encode('utf-8')
    lname = struct.pack('<B', len(read_name))
    mapq = struct.pack('<B', rec.mapping_quality)
    bin_val = struct.pack('<H', rec.bin)
    # Number of CIGAR operations as unsigned short.
    ncigar = struct.pack('<H', len(rec.cigar) * 4 if rec.cigar else 0)
    flags = struct.pack('<H', rec.flag)
    seq_len = struct.pack('<I', rec.query_length)
    next_refid = struct.pack('<i', rec.next_reference_id)
    next_pos = struct.pack('<i', rec.next_reference_start)
    template_len = struct.pack('<i', rec.template_length)

    # Append fixed fields into their column buffers.
    columns["RefID"].extend(refid)
    columns["Pos"].extend(pos)
    columns["LName"].extend(lname)
    columns["Mapq"].extend(mapq)
    columns["Bin"].extend(bin_val)
    columns["NCigar"].extend(ncigar)
    columns["Flags"].extend(flags)
    columns["SequenceLength"].extend(seq_len)
    columns["NextRefID"].extend(next_refid)
    columns["NextPos"].extend(next_pos)
    columns["TemplateLength"].extend(template_len)

    # Variable fields:
    # ReadName: store raw bytes of the read name.
    columns["ReadName"].extend(read_name)

    # RawCigar: pack each CIGAR operation as a 4-byte integer.
    if rec.cigar:
        for op, length in rec.cigar:
            # Combine length and op into a 32-bit int as in the BAM spec.
            cigar_val = (length << 4) | op
            columns["RawCigar"].extend(struct.pack('<I', cigar_val))
    # RawSequence: pack two bases per byte.
    seq = rec.query_sequence
    raw_seq = bytearray()
    base_encoding = {'=':0, 'A':1, 'C':2, 'M':3, 'G':4, 'R':5, 'S':6, 'V':7,
                     'T':8, 'W':9, 'Y':10, 'H':11, 'K':12, 'D':13, 'B':14, 'N':15}
    i = 0
    while i < len(seq):
        first = base_encoding.get(seq[i], 15)
        second = base_encoding.get(seq[i+1], 0) if i+1 < len(seq) else 0
        raw_seq.append((first << 4) | second)
        i += 2
    columns["RawSequence"].extend(raw_seq)
    # Also record the length of the raw sequence in a separate column.
    columns["RawSeqLen"].extend(struct.pack('<I', len(raw_seq)))

    # RawQual: quality scores (as bytes).
    qual = bytes(rec.query_qualities) if rec.query_qualities else b'\xff' * rec.query_length
    columns["RawQual"].extend(qual)

    # RawTags: here we simply serialize the tags as JSON (for illustration).
    # In a production version  would pack the binary representation exactly.
    tags_json = json.dumps(rec.get_tags())
    tags_bytes = tags_json.encode('utf-8')
    columns["RawTags"].extend(tags_bytes)
    columns["RawTagsLen"].extend(struct.pack('<I', len(tags_bytes)))

    # (Optionally: flush any column that grows too large.)
    for field, buf in columns.items():
        if len(buf) > BLOCK_SIZE_LIMIT:
            # Flush this column’s block to file later.
            pass  # We delay flushing until writing out the file.

# --- Write out the GBAM file ---
with open(gbam_output_path, "wb") as out_f:
    # Reserve first 1000 bytes for file info header.
    header_placeholder = b'\x00' * 1000
    out_f.write(header_placeholder)

    # Write each column’s data (flush any remaining data as one block per column).
    for field in fields_mapping:
        flush_column(field, out_f)

    # Record the current file position (meta start position).
    meta_seekpos = out_f.tell()

    # Create file meta data.
    file_meta = {
        "field_to_meta": field_meta,  # mapping of field names to their block meta info
        "sam_header": list(sam_header_bytes),  # store as a list of integers (or base64 encode)
        "name_to_ref_id": ref_seqs
    }
    meta_json = json.dumps(file_meta).encode('utf-8')
    out_f.write(meta_json)

    # Now update the 1000-byte header.
    file_info = {
        "magic": "GBAM",
        "gbam_version": [1, 0],
        "seekpos": meta_seekpos,
        "crc32": 0,  #  could compute CRC32 over meta_json if desired
        "is_sorted": True,
        "creation_command": "python bam_to_gbam.py"
    }
    file_info_json = json.dumps(file_info).encode('utf-8')
    header = file_info_json.ljust(1000, b'\x00')
    out_f.seek(0)
    out_f.write(header)

print("GBAM file successfully created: {gbam_output_path}")
