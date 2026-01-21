import argparse
import json
import struct
import lz4.block
import pysam
import zstandard as zstd
import brotli
import shutil
import sys
import os
import psutil
from collections import OrderedDict
from enum import Enum
from preprocess import preprocess_readnames_for_tokenizer

# --- Settings ---
BATCH_SIZE = 1_000_000
FIXED_FIELD_SIZE = {
    "RefID": 4, "Pos": 4, "LName": 1, "Mapq": 1, "Bin": 2, "NCigar": 2,
    "Flags": 2, "SequenceLength": 4, "NextRefID": 4, "NextPos": 4,
    "TemplateLength": 4, "RawSeqLen": 2, "RawTagsLen": 4
}

# --- Compression codec options ---
class CompressionCodec(str, Enum):
    LZ4 = "lz4"
    ZSTD = "zstd"
    BROTLI = "brotli"

# --- Field mapping and compression ---
column_compression = {field: CompressionCodec.BROTLI for field in [
    "RefID", "Pos", "LName", "Mapq", "Bin", "NCigar", "Flags", "SequenceLength",
    "NextRefID", "NextPos", "TemplateLength", "RawSeqLen", "RawTagsLen",
    "ReadName", "RawCigar", "RawSequence", "RawQual", "RawTags"]}

fields_mapping = OrderedDict({
    "RefID": None, "Pos": None, "LName": None, "Mapq": None, "Bin": None,
    "NCigar": None, "Flags": None, "SequenceLength": None, "NextRefID": None,
    "NextPos": None, "TemplateLength": None, "RawSeqLen": None,
    "RawTagsLen": None, "ReadName": "LName", "RawCigar": "NCigar",
    "RawSequence": "RawSeqLen", "RawQual": "SequenceLength",
    "RawTags": "RawTagsLen"
})

# --- Output state ---
columns = {field: bytearray() for field in fields_mapping}
field_meta = {field: [] for field in fields_mapping}
readnames_list = []
record_count = 0

# --- Argument parsing ---
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", required=True)
parser.add_argument("-o", "--output", required=True)
parser.add_argument("--compression-config")
args = parser.parse_args()

if args.compression_config:
    with open(args.compression_config) as f:
        raw_map = json.load(f)
        column_compression = {k: CompressionCodec(raw_map[k]) for k in raw_map}

bam_file = pysam.AlignmentFile(args.input, "rb")
sam_header_bytes = bam_file.text.encode("utf-8")
ref_seqs = [(sq['SN'], sq['LN']) for sq in bam_file.header.get('SQ', [])]
out_f = open(args.output, "wb")
out_f.write(b"\x00" * 1000)
gbam_dir = os.path.dirname(args.output) or "."

# --- Flush helper ---
def flush_column(field, out_f, records_in_block):
    buf = columns[field]
    if not buf:
        return
    uncompressed_data = bytes(buf)
    codec = column_compression.get(field, CompressionCodec.BROTLI)

    if codec == CompressionCodec.LZ4:
        compressed_data = lz4.block.compress(uncompressed_data, store_size=False)
    elif codec == CompressionCodec.ZSTD:
        compressed_data = zstd.ZstdCompressor(level=3).compress(uncompressed_data)
    elif codec == CompressionCodec.BROTLI:
        compressed_data = brotli.compress(uncompressed_data, quality=11)
    else:
        raise ValueError(f"Unknown codec {codec}")

    seekpos = out_f.tell()
    out_f.write(compressed_data)

    meta = {
        "seekpos": seekpos,
        "numitems": records_in_block,
        "block_size": len(compressed_data),
        "uncompressed_size": len(uncompressed_data),
        "codec": codec
    }
    field_meta[field].append(meta)
    columns[field].clear()

# --- Process batch ---
def process_batch(batch):
    global record_count
    batch.sort(key=lambda r: (r.reference_id, r.mapping_quality, r.flag, r.reference_start, r.query_name))
    for rec in batch:
        qname = rec.query_name
        readnames_list.append(qname)
        encoded_name = qname.encode('utf-8')

        columns["RefID"].extend(struct.pack('<i', rec.reference_id))
        columns["Pos"].extend(struct.pack('<i', rec.reference_start))
        columns["LName"].extend(struct.pack('<B', len(encoded_name)))
        columns["Mapq"].extend(struct.pack('<B', rec.mapping_quality))
        columns["Bin"].extend(struct.pack('<H', rec.bin))
        columns["NCigar"].extend(struct.pack('<H', len(rec.cigar) * 4 if rec.cigar else 0))
        columns["Flags"].extend(struct.pack('<H', rec.flag))
        columns["SequenceLength"].extend(struct.pack('<I', rec.query_length))
        columns["NextRefID"].extend(struct.pack('<i', rec.next_reference_id))
        columns["NextPos"].extend(struct.pack('<i', rec.next_reference_start))
        columns["TemplateLength"].extend(struct.pack('<i', rec.template_length))
        columns["ReadName"].extend(encoded_name)

        if rec.cigar:
            for op, length in rec.cigar:
                columns["RawCigar"].extend(struct.pack('<I', (length << 4) | op))

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
        columns["RawSeqLen"].extend(struct.pack('<H', len(raw_seq)))

        qual = bytes(rec.query_qualities) if rec.query_qualities else b'\xff' * rec.query_length
        columns["RawQual"].extend(qual)

        tags_bytes = json.dumps(rec.get_tags()).encode('utf-8')
        columns["RawTags"].extend(tags_bytes)
        columns["RawTagsLen"].extend(struct.pack('<I', len(tags_bytes)))

        record_count += 1

    for field in fields_mapping:
        flush_column(field, out_f, len(batch))

# --- Main loop ---
batch = []
for rec in bam_file.fetch(until_eof=True):
    batch.append(rec)
    if len(batch) >= BATCH_SIZE:
        process_batch(batch)
        batch.clear()
        mem = psutil.Process(os.getpid()).memory_info().rss / (1024 ** 2)
        print(f"[INFO] Processed {record_count} records. Memory usage: {mem:.2f} MB")

if batch:
    process_batch(batch)

compressed_readnames, readname_vocab = preprocess_readnames_for_tokenizer(columns, readnames_list)
columns["ReadName"] = bytearray(compressed_readnames)

meta_json = json.dumps({
    "field_to_meta": field_meta,
    "sam_header": list(sam_header_bytes),
    "name_to_ref_id": ref_seqs,
    "readname_vocab": readname_vocab
}).encode('utf-8')
compressed_meta = brotli.compress(meta_json, quality=11)
meta_seekpos = out_f.tell()
out_f.write(compressed_meta)

file_info = {
    "magic": "GBAM",
    "gbam_version": [1, 0],
    "seekpos": meta_seekpos,
    "crc32": 0,
    "is_sorted": True,
    "creation_command": "python bam_to_gbam.py"
}
out_f.seek(0)
out_f.write(json.dumps(file_info).encode('utf-8').ljust(1000, b'\x00'))
out_f.close()
print(f"[INFO] GBAM file written to {args.output}")
