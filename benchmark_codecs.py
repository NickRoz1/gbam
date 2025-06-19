import json
import subprocess
import os
import csv

# Config
CODEC_MAP_PATH = "codec_map.json"
INPUT_FILE = "test_data/little.bam"
OUTPUT_FILE = "test_data/compression_test.gbam"
CSV_OUTPUT = "codec_benchmark_results.csv"
CODECS = ["Zstd", "Gzip", "Brotli", "Lz4"]

def load_codec_map():
    with open(CODEC_MAP_PATH, "r") as f:
        return json.load(f)

def save_codec_map(codec_map):
    with open(CODEC_MAP_PATH, "w") as f:
        json.dump(codec_map, f, indent=2)

def run_encoder():
    result = subprocess.run([
        "./target/release/gbam_binary",
        "-c", INPUT_FILE,
        "-o", OUTPUT_FILE
    ], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode != 0:
        print("Binary execution failed:")
        print(result.stderr.decode())
        raise RuntimeError("Binary execution failed.")

def measure_size(filepath):
    return os.path.getsize(filepath)

def main():
    base_map = load_codec_map()
    fields = list(base_map.keys())

    # Step 1: Baseline with all NoCompression
    baseline_map = {f: "NoCompression" for f in base_map}
    save_codec_map(baseline_map)
    print("📏 Measuring baseline (all NoCompression)...")
    run_encoder()
    baseline_size = measure_size(OUTPUT_FILE)
    print(f"📦 Baseline GBAM size: {baseline_size} bytes\n")

    # Step 2: Test each field with each codec
    results = {}  # Dict: field -> {codec -> %saved}

    for field in fields:
        results[field] = {}

        for codec in CODECS:
            # Set only this field to codec, others to NoCompression
            test_map = {f: "NoCompression" for f in base_map}
            test_map[field] = codec
            save_codec_map(test_map)

            print(f"Testing {field} with {codec}...")
            run_encoder()

            compressed_size = measure_size(OUTPUT_FILE)
            space_saved = baseline_size - compressed_size
            space_saved_percent = (space_saved / baseline_size) * 100

            results[field][codec] = f"{space_saved_percent:.2f}"

        # Optional: Reset field (not strictly needed due to overwrite above)
        base_map[field] = "NoCompression"
        save_codec_map(base_map)

    # Step 3: Write results to CSV
    with open(CSV_OUTPUT, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        header = ["Field"] + [f"{codec} %" for codec in CODECS]
        writer.writerow(header)

        for field in fields:
            row = [field] + [results[field].get(codec, "0.00") for codec in CODECS]
            writer.writerow(row)

    print(f"\nFinal results saved to: {CSV_OUTPUT}")

if __name__ == "__main__":
    main()
