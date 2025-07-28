import os
import time
import subprocess
import tempfile
import json
from pathlib import Path

# Paths to scripts
BAM_TO_GBAM_SCRIPT = "bam_to_gbam.py"
GBAM_TO_BAM_SCRIPT = "gbam_to_bam.py"

# Input/output files
ORIGINAL_BAM = "/Users/hasitha/Documents/biology/gbam/gbam_python/test_data/HG001.bam"
GBAM_OUTPUT = "test_data/output.gbam"
RECONSTRUCTED_BAM = "test_data/reconstructed.bam"

# Compression test cases
test_cases = {
    "all_lz4": {field: "lz4" for field in [
        "RefID", "Pos", "LName", "Mapq", "Bin", "NCigar", "Flags", "SequenceLength",
        "NextRefID", "NextPos", "TemplateLength", "RawSeqLen", "RawTagsLen",
        "ReadName", "RawCigar", "RawSequence", "RawQual", "RawTags"
    ]},
    "all_zstd": {field: "zstd" for field in [
        "RefID", "Pos", "LName", "Mapq", "Bin", "NCigar", "Flags", "SequenceLength",
        "NextRefID", "NextPos", "TemplateLength", "RawSeqLen", "RawTagsLen",
        "ReadName", "RawCigar", "RawSequence", "RawQual", "RawTags"
    ]},
    # "all_brotli": {field: "brotli" for field in [
    #     "RefID", "Pos", "LName", "Mapq", "Bin", "NCigar", "Flags", "SequenceLength",
    #     "NextRefID", "NextPos", "TemplateLength", "RawSeqLen", "RawTagsLen",
    #     "ReadName", "RawCigar", "RawSequence", "RawQual", "RawTags"
    # ]},
    "zstd_fixed_brotli_var": {
        field: "zstd" if index is None else "brotli"
        for field, index in [
            ("RefID", None), ("Pos", None), ("LName", None), ("Mapq", None), ("Bin", None),
            ("NCigar", None), ("Flags", None), ("SequenceLength", None), ("NextRefID", None),
            ("NextPos", None), ("TemplateLength", None), ("RawSeqLen", None), ("RawTagsLen", None),
            ("ReadName", "LName"), ("RawCigar", "NCigar"), ("RawSequence", "RawSeqLen"),
            ("RawQual", "SequenceLength"), ("RawTags", "RawTagsLen")
        ]
    },
    "brotli_fixed_zstd_var": {
        field: "brotli" if index is None else "zstd"
        for field, index in [
            ("RefID", None), ("Pos", None), ("LName", None), ("Mapq", None), ("Bin", None),
            ("NCigar", None), ("Flags", None), ("SequenceLength", None), ("NextRefID", None),
            ("NextPos", None), ("TemplateLength", None), ("RawSeqLen", None), ("RawTagsLen", None),
            ("ReadName", "LName"), ("RawCigar", "NCigar"), ("RawSequence", "RawSeqLen"),
            ("RawQual", "SequenceLength"), ("RawTags", "RawTagsLen")
        ]
    },
    "hybrid": {
        "RefID": "lz4", "Pos": "lz4", "LName": "lz4", "Mapq": "lz4", "Bin": "lz4",
        "NCigar": "lz4", "Flags": "lz4", "SequenceLength": "lz4", "NextRefID": "lz4",
        "NextPos": "zstd", "TemplateLength": "zstd", "RawSeqLen": "zstd", "RawTagsLen": "zstd",
        "ReadName": "brotli", "RawCigar": "zstd", "RawSequence": "zstd", "RawQual": "zstd",
        "RawTags": "brotli"
    }
}

def run_command(cmd):
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(result.stderr)
        raise RuntimeError(f"Command failed: {' '.join(cmd)}")
    return result.stdout

def benchmark_case(name, compression_map):
    print(f"\n== Running test: {name} ==")
    bam_size = os.path.getsize(ORIGINAL_BAM)

    # Write compression config
    with tempfile.NamedTemporaryFile(mode='w', suffix=".json", delete=False) as tmp:
        json.dump(compression_map, tmp)
        tmp.flush()
        config_path = tmp.name

    # 1. BAM → GBAM
    run_command(["python3", BAM_TO_GBAM_SCRIPT, "-i", ORIGINAL_BAM, "-o", GBAM_OUTPUT, "--compression-config", config_path])
    gbam_size = os.path.getsize(GBAM_OUTPUT)

    # 2. GBAM → BAM (measure time)
    start = time.time()
    run_command(["python3", GBAM_TO_BAM_SCRIPT, "-i", GBAM_OUTPUT, "-o", RECONSTRUCTED_BAM])
    duration = time.time() - start
    throughput = gbam_size / (1024 * 1024) / duration

    os.remove(config_path)

    return {
        "test": name,
        "bam_size_mb": round(bam_size / (1024 * 1024), 2),
        "gbam_size_mb": round(gbam_size / (1024 * 1024), 2),
        "decompression_time_s": round(duration, 2),
        "throughput_mb_s": round(throughput, 2)
    }

def main():
    results = []
    for name, compression_map in test_cases.items():
        result = benchmark_case(name, compression_map)
        results.append(result)

    print("\n=== Results ===")
    for r in results:
        print(json.dumps(r, indent=2))

    # Optional: write CSV
    with open("test_data/compression_benchmark.csv", "w") as f:
        f.write("test,bam_size_mb,gbam_size_mb,decompression_time_s,throughput_mb_s\n")
        for r in results:
            f.write(f"{r['test']},{r['bam_size_mb']},{r['gbam_size_mb']},{r['decompression_time_s']},{r['throughput_mb_s']}\n")

if __name__ == "__main__":
    main()
