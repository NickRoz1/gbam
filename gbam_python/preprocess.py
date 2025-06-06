# --- Delta + RLE Encoding for RawQual ---
from itertools import groupby
import struct
import re
from typing import List, Tuple
import brotli
import json

# Delta encode RawQual

def delta_encode_rawqual(qual_bytes):
    if not qual_bytes:
        return b''
    deltas = [qual_bytes[0]]
    for i in range(1, len(qual_bytes)):
        delta = (qual_bytes[i] - qual_bytes[i - 1]) % 256  # wrap-around delta
        deltas.append(delta)
    return bytes(deltas)

# Run-length encode RawQual

def rle_encode_rawqual(qual_bytes):
    rle_encoded = bytearray()
    for val, group in groupby(qual_bytes):
        run = list(group)
        count = len(run)
        while count > 255:
            rle_encoded.append(255)
            rle_encoded.append(val)
            count -= 255
        rle_encoded.append(count)
        rle_encoded.append(val)
    return bytes(rle_encoded)

# --- Usage Example ---
def preprocess_rawqual_for_brotli(columns):
    raw_qual = bytes(columns["RawQual"])

    # Option 1: Delta encoding
    # encoded = delta_encode_rawqual(raw_qual)

    # Option 2: Run-length encoding (more effective for repetitive values)
    encoded = rle_encode_rawqual(raw_qual)

    print(f"[INFO] Preprocessed RawQual from {len(raw_qual)} to {len(encoded)} bytes using RLE")
    columns["RawQual"] = bytearray(encoded)

# --- Tokenization + Delta Encoding for ReadName ---
def tokenize_readnames(readname_list):
    token_columns = []
    max_tokens = max(name.count(":") + 1 for name in readname_list)
    for _ in range(max_tokens):
        token_columns.append([])

    for name in readname_list:
        tokens = name.lstrip("@").split(":")
        for i in range(len(tokens)):
            token_columns[i].append(tokens[i])

    return token_columns

def delta_encode_column(col):
    if not col:
        return []
    try:
        col = [int(x) for x in col]
        deltas = [col[0]]
        for i in range(1, len(col)):
            deltas.append(col[i] - col[i - 1])
        return deltas
    except ValueError:
        return col  # non-numeric, skip delta

def preprocess_readnames_for_brotli(columns):
    raw_names = columns["ReadName"].decode("utf-8").split("\n")
    token_cols = tokenize_readnames(raw_names)
    result = []
    for col in token_cols:
        encoded = delta_encode_column(col)
        result.append(encoded)

    # Optional: save or convert to compact binary form before compression
    flat = json.dumps(result).encode("utf-8")
    print(f"[INFO] Tokenized + delta-encoded ReadName from {len(raw_names)} to {len(flat)} bytes")
    columns["ReadName"] = bytearray(flat)
    
def compress_vocab_brotli(vocab: List[str]) -> bytes:
    """
    Compress the vocabulary list using Brotli.
    Format: JSON string of vocab, then Brotli-compressed.
    """
    vocab_json = json.dumps(vocab).encode('utf-8')
    compressed_vocab = brotli.compress(vocab_json, quality=11)
    print(f"[INFO] Compressed vocab size: {len(compressed_vocab)} bytes from {len(vocab_json)} raw bytes")
    return compressed_vocab

def preprocess_readnames_for_tokenizer(columns: dict, input_read_names: List) -> Tuple[bytes, List[str]]:
    """
    Tokenizes and compresses read names using ':' as the sole separator.
    Uses binary serialization and Brotli compression.
    Returns compressed data and vocabulary list for metadata storage.
    """

    def tokenize_read_names(read_names: List[str]) -> Tuple[List[List[str]], List[str]]:
        tokens = []
        all_parts = []

        for name in read_names:
            parts = name.split(':')  # Only split on ':'
            tokens.append(parts)
            all_parts.extend(parts)

        # Add __PAD__ token for safe decoding
        unique_tokens = set(all_parts)
        unique_tokens.add("__PAD__")
        vocab = list(sorted(unique_tokens))
        return tokens, vocab

    def build_token_indices(tokens: List[List[str]], vocab: List[str]) -> Tuple[List[List[int]], int]:
        token_to_index = {token: idx for idx, token in enumerate(vocab)}
        pad_id = token_to_index["__PAD__"]
        indexed_tokens = [[token_to_index[token] for token in row] for row in tokens]
        return indexed_tokens, pad_id

    def compress_tokenized_names_binary(indexed_tokens: List[List[int]], pad_id: int) -> bytes:
        output = bytearray()
        num_records = len(indexed_tokens)
        max_len = max(len(toks) for toks in indexed_tokens)
        padded = [toks + [pad_id] * (max_len - len(toks)) for toks in indexed_tokens]

        output.extend(struct.pack('<I', num_records))
        output.extend(struct.pack('<I', max_len))

        for row in padded:
            for token_id in row:
                output.extend(struct.pack('<H', token_id))  # uint16

        return bytes(output)

    # Extract read names
    # read_names = [rec.query_name for rec in all_records]
    tokens, vocab = tokenize_read_names(input_read_names)
    indexed_tokens, pad_id = build_token_indices(tokens, vocab)
    binary_serialized = compress_tokenized_names_binary(indexed_tokens, pad_id)
    compressed_bytes = brotli.compress(binary_serialized, quality=11)

    print(f"[INFO] Tokenizer + Brotli compressed ReadName column from {sum(len(name.encode()) for name in input_read_names)} to {len(compressed_bytes)} bytes.")
    print(f"[INFO] ReadName tokenizer vocab: {len(vocab)} tokens, {sum(len(t.encode()) for t in vocab)} bytes total")

    columns["ReadName"] = bytearray(compressed_bytes)
    return compressed_bytes, vocab

def decompress_readnames_from_tokenizer(compressed_data: bytes, vocab: List[str]) -> List[str]:
    """
    Decompresses Brotli-compressed tokenized read names and reconstructs original names.
    Assumes __PAD__ token was used during encoding for row padding.
    """
    # Step 1: Brotli decompress
    binary_data = brotli.decompress(compressed_data)

    # Step 2: Read header values
    num_records = struct.unpack_from('<I', binary_data, 0)[0]
    max_len = struct.unpack_from('<I', binary_data, 4)[0]

    # Step 3: Extract all token indices
    offset = 8
    total_tokens = num_records * max_len
    token_ids = struct.unpack_from(f'<{total_tokens}H', binary_data, offset)

    # Step 4: Reconstruct token rows
    indexed_tokens = [list(token_ids[i * max_len:(i + 1) * max_len]) for i in range(num_records)]

    # Step 5: Find the PAD token index
    try:
        pad_token_index = vocab.index("__PAD__")
    except ValueError:
        pad_token_index = None  # Fallback if __PAD__ is missing

    # Step 6: Decode token rows back to strings
    read_names = []
    for row in indexed_tokens:
        if pad_token_index is not None:
            tokens = [vocab[idx] for idx in row if idx != pad_token_index]
        else:
            tokens = [vocab[idx] for idx in row]
        read_names.append(":".join(tokens))

    return read_names

def encode_sequence_2bit(seq: str) -> bytes:
    """Encode a nucleotide sequence into 2-bit packed bytes."""
    base2bit = {
        'A': 0b00,
        'C': 0b01,
        'G': 0b10,
        'T': 0b11,
    }
    encoded = 0
    bits_filled = 0
    buffer = bytearray()

    for base in seq.upper():
        val = base2bit.get(base, 0b00)  # default to 'A' if unknown
        encoded = (encoded << 2) | val
        bits_filled += 2

        if bits_filled == 8:
            buffer.append(encoded)
            encoded = 0
            bits_filled = 0

    if bits_filled > 0:
        encoded <<= (8 - bits_filled)  # pad remaining bits
        buffer.append(encoded)

    return bytes(buffer)

def preprocess_rawqual_fqz_style(raw_quals, read_lengths):
    """
    Emits (delta, pos_bin) as two separate bytes per base.
    """
    output = bytearray()
    i = 0
    for length in read_lengths:
        prev = 0
        for pos in range(length):
            q = raw_quals[i]
            delta = (q - prev) & 0xFF  # full byte delta
            pos_bin = min(pos // 4, 255)  # safe cap
            output.append(delta)
            output.append(pos_bin)
            prev = q
            i += 1
    return output

