#!/bin/bash
set -e  # Exit on error

# Create unique temp files in /dev/shm
t1=$(mktemp /dev/shm/gbam_output.XXXXXX)
t2=$(mktemp /dev/shm/tsamtools_output.XXXXXX)

# Ensure cleanup on script exit or crash
trap 'rm -f "$t1" "$t2"' EXIT

./main output.gbam < ../test_data/100recs.bam
./main output.gbam > "$t1"
samtools view ../test_data/100recs.bam > "$t2"

if diff -q "$t1" "$t2" > /dev/null; then
    exit 0
else
    echo "Test failed: Outputs differ!"
    exit 1
fi