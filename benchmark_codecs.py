import json
import subprocess
import os
import csv
import argparse

# Constants
CODEC_MAP_PATH = "codec_map.json"
CSV_OUTPUT = "codec_benchmark_results.csv"
CODECS = ["Zstd", "Gzip", "Brotli", "Lz4"]

def load_codec_map():
    with open(CODEC_MAP_PATH, "r") as f:
        return json.load(f)

def save_codec_map(codec_map):
    with open(CODEC_MAP_PATH, "w") as f:
        json.dump(codec_map, f, indent=2)

def run_encoder(input_file, output_file):
    result = subprocess.run([
        "./target/release/gbam_binary",
        "-c", input_file,
        "-o", output_file
    ], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode != 0:
        print("Binary execution failed:")
        print(result.stderr.decode())
        raise RuntimeError("Binary execution failed.")

def measure_size(filepath):
    return os.path.getsize(filepath)

def main():
    parser = argparse.ArgumentParser(description="Benchmark per-field compression")
    parser.add_argument("-i", "--input", required=True, help="Input BAM file path")
    parser.add_argument("-o", "--output", required=True, help="Output GBAM file path")
    args = parser.parse_args()

    input_file = args.input
    output_file = args.output

    base_map = load_codec_map()
    fields = list(base_map.keys())

    # Step 1: Baseline
    baseline_map = {f: "NoCompression" for f in base_map}
    save_codec_map(baseline_map)
    print("Measuring baseline (all NoCompression)...")
    run_encoder(input_file, output_file)
    baseline_size = measure_size(output_file)
    print(f"Baseline GBAM size: {baseline_size} bytes\n")

    # Step 2: Test each field
    results = {}  # field -> {codec -> %saved}
    for field in fields:
        results[field] = {}
        for codec in CODECS:
            test_map = {f: "NoCompression" for f in base_map}
            test_map[field] = codec
            save_codec_map(test_map)

            print(f"Testing {field} with {codec}...")
            run_encoder(input_file, output_file)

            compressed_size = measure_size(output_file)
            space_saved = baseline_size - compressed_size
            space_saved_percent = (space_saved / baseline_size) * 100

            results[field][codec] = f"{space_saved_percent:.2f}"

    # Step 3: Write CSV
    with open(CSV_OUTPUT, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        header = ["Field"] + [f"{codec} %" for codec in CODECS]
        writer.writerow(header)

        for field in fields:
            row = [field] + [results[field].get(codec, "0.00") for codec in CODECS]
            writer.writerow(row)

    print(f"Benchmark complete. Results saved to: {CSV_OUTPUT}")

if __name__ == "__main__":
    main()
