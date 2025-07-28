# GBAM Compression Benchmark

This prototype evaluates different compression strategies for converting BAM to GBAM files, analyzing:

- Original BAM size
- Compressed GBAM size
- Decompression time
- Throughput (MB/s during GBAM → BAM conversion)

## How to Run the Benchmark

### 1. **Prerequisites**
Ensure the following are installed:
- Python 3
- `pysam`, `lz4`, `zstandard`, `brotli`

Install with:
```bash
pip install pysam lz4 zstandard brotli
```

Run the test:
```bash
python benchmark_test.py
```

## Benchmark Summary

| Test Name              | BAM Size (MB) | GBAM Size (MB) | Decompression Time (s) | Throughput (MB/s)  |
|------------------------|---------------|----------------|------------------------|--------------------|
| all_lz4                | 49.33         | 68.37          | 54.08                  | 1.26               |
| all_zstd               | 49.33         | 41.35          | 52.01                  | 0.79               |
| all_brotli             | 49.33         | 36.64          | 53.16                  | 0.69               |
| zstd_fixed_brotli_var  | 49.33         | 37.02          | 53.10                  | 0.70               |
| brotli_fixed_zstd_var  | 49.33         | 40.97          | 52.88                  | 0.77               |
| hybrid                 | 49.33         | 39.40          | 52.38                  | 0.75               |

## Test Descriptions

- **all_lz4**: All fields compressed using LZ4 (fastest decompression, largest file).
- **all_zstd**: All fields compressed using Zstandard (balanced).
- **all_brotli**: All fields compressed using Brotli (smallest file, slowest).
- **zstd_fixed_brotli_var**: Fixed-width numeric fields with Zstd, text-based fields with Brotli.
- **brotli_fixed_zstd_var**: Opposite of the above.
- **hybrid**: Hand-tuned mix for best balance.
