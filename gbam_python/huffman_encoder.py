# --- Huffman Utilities ---
from collections import Counter, namedtuple
import heapq

class HuffmanNode(namedtuple("Node", ["char", "freq", "left", "right"])):
    def __lt__(self, other):
        return self.freq < other.freq

def build_huffman_tree(data):
    freq_map = Counter(data)
    heap = [HuffmanNode(char=c, freq=f, left=None, right=None) for c, f in freq_map.items()]
    heapq.heapify(heap)
    while len(heap) > 1:
        a = heapq.heappop(heap)
        b = heapq.heappop(heap)
        heapq.heappush(heap, HuffmanNode(None, a.freq + b.freq, a, b))
    return heap[0] if heap else None

def build_codes(node, prefix='', codebook=None):
    if codebook is None:
        codebook = {}
    if node.char is not None:
        codebook[node.char] = prefix
    else:
        build_codes(node.left, prefix + '0', codebook)
        build_codes(node.right, prefix + '1', codebook)
    return codebook

def huffman_encode(data, codebook):
    return ''.join(codebook[byte] for byte in data)

def bitstring_to_bytes(bitstring):
    pad_len = (8 - len(bitstring) % 8) % 8
    bitstring += '0' * pad_len
    return bytes(int(bitstring[i:i+8], 2) for i in range(0, len(bitstring), 8)), pad_len

# --- Huffman Compression for RawQual ---
def compress_rawqual_with_huffman(columns, field_meta):
    raw_qual_bytes = bytes(columns["RawQual"])
    huff_tree = build_huffman_tree(raw_qual_bytes)

    if huff_tree:
        huff_codes = build_codes(huff_tree)
        encoded_bits = huffman_encode(raw_qual_bytes, huff_codes)
        compressed_rawqual_bytes, pad_len = bitstring_to_bytes(encoded_bits)

        print(f"[INFO] Huffman compressed RawQual from {len(raw_qual_bytes)} to {len(compressed_rawqual_bytes)} bytes")

        # Replace RawQual column with compressed data
        columns["RawQual"] = bytearray(compressed_rawqual_bytes)

        # Add Huffman metadata for later decoding (optional)
        field_meta["RawQual"].append({
            "huffman_codes": {k: v for k, v in huff_codes.items()},
            "padding": pad_len,
            "original_size": len(raw_qual_bytes),
            "compressed_size": len(compressed_rawqual_bytes)
        })
    else:
        print("[WARN] RawQual is empty or has only one unique value — Huffman skipped.")
