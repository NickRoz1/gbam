#! /usr/bin/env python3

import argparse
import subprocess
import os
import sys
from pathlib import Path
import tempfile
import time

timer = 'echo 3 > /proc/sys/vm/drop_caches ; /usr/bin/time -v'

samtools_bam_results = {}
samtools_cram_results = {}
sambamba_results = {}
gbam_results = {}
gfa_inject_results = {}

keywords = ["User", "System", "Percent", "Elapsed", "Maximum"]
newline = '\n'
# my_format = lambda name, out: f"## {name}\n ```\n{out.decode('utf-8')}```\n\n"
my_format = lambda name, out: f"## {name}\n ```\n{newline.join([s.strip() for s in out.decode('utf-8').splitlines() if any(word in s for word in keywords)])}\n```\n\n"
total_time_format = lambda name, out: f"## {name}\n ```\n{out}s\n```\n\n"

def gbam(bin_path, bam_path, res_path):
    start = time.time()

    gbam_path = res_path/"index_sort_out.gbam"
    depth_gbam_out = res_path/"gbam_depth_out.txt"

    gbam_results["index-sort"]  = subprocess.check_output([f"{timer} {bin_path} {bam_path} -c -s --sort-temp-mode ram -o {gbam_path} --index-sort"], shell=True, stderr=subprocess.STDOUT)
    gbam_results["flagstat"]  = subprocess.check_output([f"{timer} {bin_path} --flagstat {gbam_path}"], shell=True, stderr=subprocess.STDOUT)
    gbam_results["depth"] = subprocess.check_output([f"{timer} {bin_path} --depth {gbam_path} --thread-num 12 --index-file {gbam_path.with_suffix('.gbam.gbai')} > {depth_gbam_out}"], shell=True, stderr=subprocess.STDOUT)

    print(f"Completed GBAM benchmarking, took: {time.time()-start} seconds")

def samtools_bam(bin_path, bam_path, res_path):
    start = time.time()

    samtools_path = res_path/"sort_out.samtools.bam"
    depth_samtools_out = res_path/"samtools_depth_out.txt"

    samtools_bam_results["sort"] = subprocess.check_output([f"{timer} {bin_path} sort -@ 8 {bam_path} -o {samtools_path}"], shell=True, stderr=subprocess.STDOUT)
    samtools_bam_results["flagstat"] = subprocess.check_output([f"{timer} {bin_path} flagstat -@ 8 {samtools_path}"], shell=True, stderr=subprocess.STDOUT)
    samtools_bam_results["depth"] = subprocess.check_output([f"{timer} {bin_path} depth -@ 8 {samtools_path} > {depth_samtools_out}"], shell=True, stderr=subprocess.STDOUT)

    print(f"Completed SAMTOOLS BAM benchmarking, took: {time.time()-start} seconds")

def samtools_cram(bin_path, bam_path, res_path):
    print("Creating cram file from bam file.")
    test_cram_file = tempfile.NamedTemporaryFile(suffix=".cram")
    subprocess.check_output([f"{bin_path} view -C -o {test_cram_file.name} {bam_path}"], shell=True, stderr=subprocess.STDOUT)
    print("Starting benchmark using created cram file.")

    start = time.time()
    
    samtools_path = res_path/"sort_out.samtools.cram"
    depth_samtools_out = res_path/"samtools_depth_out.txt"

    samtools_cram_results["sort"] = subprocess.check_output([f"{timer} {bin_path} sort -@ 8 {test_cram_file.name} -o {samtools_path}"], shell=True, stderr=subprocess.STDOUT)
    samtools_cram_results["flagstat"] = subprocess.check_output([f"{timer} {bin_path} flagstat -@ 8 {samtools_path}"], shell=True, stderr=subprocess.STDOUT)
    samtools_cram_results["depth"] = subprocess.check_output([f"{timer} {bin_path} depth -@ 8 {samtools_path} > {depth_samtools_out}"], shell=True, stderr=subprocess.STDOUT)

    print(f"Completed SAMTOOLS CRAM benchmarking, took: {time.time()-start} seconds")

def sambamba(bin_path, disable_sambamba_depth, bam_path, res_path):
    start = time.time()

    sambamba_path = res_path/"sort_out.sambamba.bam"
    depth_sambamba_out = res_path/"sambamba_depth_out.txt"

    sambamba_results["sort"] = subprocess.check_output([f"{timer} {bin_path} sort -t 8 -m 90GB {bam_path} -o {sambamba_path}"], shell=True, stderr=subprocess.STDOUT)
    sambamba_results["flagstat"] = subprocess.check_output([f"{timer} {bin_path} flagstat -t 8 {sambamba_path}"], shell=True, stderr=subprocess.STDOUT)
    if not disable_sambamba_depth:
        sambamba_results["depth"] = subprocess.check_output([f"{timer} {bin_path} depth base -t 8 {sambamba_path} > {depth_sambamba_out}"], shell=True, stderr=subprocess.STDOUT)
    else:
        sambamba_results["depth"] = subprocess.check_output([f"{timer} echo 'Test is disabled'"], shell=True, stderr=subprocess.STDOUT)
    print(f"Completed SAMBAMBA benchmarking, took: {time.time()-start} seconds")

def gfainject(bin_path, gbam_path, gfa_path, bam_path, res_path):
    gfainject_from_bam_path = res_path/"injection_from_bam.gfainject.gaf"
    start = time.time()
    gfa_inject_results["bam_inject"] = subprocess.check_output([f"{timer} {bin_path} --gfa {gfa_path} --bam {bam_path }> {gfainject_from_bam_path}"], shell=True, stderr=subprocess.STDOUT)
    print(f"Completed injection from BAM, took: {time.time()-start} seconds")

    gfainject_from_gbam_path = res_path/"injection_from_gbam.gfainject.gaf"
    start = time.time()
    gfa_inject_results["gbam_inject"] = subprocess.check_output([f"{timer} {bin_path} --gfa {gfa_path} --gbam {gbam_path }> {gfainject_from_gbam_path}"], shell=True, stderr=subprocess.STDOUT)
    print(f"Completed injection from GBAM, took: {time.time()-start} seconds")

# python3 tests/benchmark.py --gbam_bin target/release/gbam_binary --bam_file test_data/little.bam --result_dir benchmarking --samtools_bin /usr/local/bin/samtools --sambamba_bin /usr/local/bin/sambamba-0.8.2-linux-amd64-static
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Benchmark SAMTOOLS/SAMBAMBA/GBAM")
    parser.add_argument("--samtools_bin",   help="Path to samtools binary", required=False)
    parser.add_argument("--sambamba_bin",   help="Path to sambamba binary", required=False)
    parser.add_argument("--gfainject_bin",  help="Path to gfainject binary", required=False)
    parser.add_argument("--gbam_bin",       help="Path to gbam binary", required=True)

    parser.add_argument("--bam_file",   help="Path to BAM file to perform test on.", required=True)
    parser.add_argument("--gbam_file",  help="Path to GBAM file to perform test on.", required=False)
    parser.add_argument("--gfa_file",   help="Path to GFA file to perform test on.", required=False)
    parser.add_argument("--result_dir", help="Path to the directory where to save the resulting files", required=True)
    parser.add_argument("--disable_sambamba_depth", help="Sambamba depth can be too slow", dest="disable_sambamba_depth", action='store_true')

    args = parser.parse_args()

    for p in [args.samtools_bin, args.sambamba_bin, args.gbam_bin, args.bam_file]:
        if not os.path.exists(p):
            print(f"Path {p} does not exist, terminating.")
            exit()

    if not os.path.exists(args.result_dir):
        os.makedirs(args.result_dir)
        print(f"Directory '{args.result_dir}' created.")

    with open(Path(args.result_dir)/"benchmark_result.md", 'w') as markdown_result_file:
        markdown_result_file.write("# GBAM Benchmarks\n\nGenerated with:\n")

        markdown_result_file.write("```sh\n{0}\n```\n\n".format(" ".join([x for x in sys.argv])))

        same_params = (args.bam_file, Path(args.result_dir))
        gbam(args.gbam_bin, *same_params)
        samtools_bam(args.samtools_bin, *same_params)
        samtools_cram(args.samtools_bin, *same_params)
        sambamba(args.sambamba_bin, args.disable_sambamba_depth, *same_params)
        if args.gfa_file is not None:
            for p in [args.gfainject_bin, args.gbam_file, args.gfa_file]:
                if not os.path.exists(p):
                    print(f"Path {p} does not exist, terminating.")
                    exit()

            gfainject(args.gfainject_bin, args.gbam_file, args.gfa_file, *same_params)
        
        # Print resulting file
        markdown_result_file.write("# SORTING\n\n\n")
        markdown_result_file.write(my_format("GBAM index-sort:", gbam_results["index-sort"]))
        markdown_result_file.write(my_format("SAMBAMBA sort:", sambamba_results["sort"]))
        markdown_result_file.write(my_format("SAMTOOLS BAM sort:", samtools_bam_results["sort"]))
        markdown_result_file.write(my_format("SAMTOOLS CRAM sort:", samtools_cram_results["sort"]))

        markdown_result_file.write("# FLAGSTAT\n\n\n")
        markdown_result_file.write(my_format("GBAM flagstat:", gbam_results["flagstat"]))
        markdown_result_file.write(my_format("SAMBAMBA flagstat:", sambamba_results["flagstat"]))
        markdown_result_file.write(my_format("SAMTOOLS BAM flagstat:", samtools_bam_results["flagstat"]))
        markdown_result_file.write(my_format("SAMTOOLS CRAM flagstat:", samtools_cram_results["flagstat"]))

        markdown_result_file.write("# DEPTH\n\n\n")
        markdown_result_file.write(my_format("GBAM depth:", gbam_results["depth"]))
        markdown_result_file.write(my_format("SAMBAMBA depth:", sambamba_results["depth"]))
        markdown_result_file.write(my_format("SAMTOOLS BAM depth:", samtools_bam_results["depth"]))
        markdown_result_file.write(my_format("SAMTOOLS CRAM depth:", samtools_cram_results["depth"]))

        if args.gfa_file is not None:
            markdown_result_file.write("# GFAINJECT\n\n\n")
            markdown_result_file.write(my_format("From GBAM:", gfa_inject_results["gbam_inject"]))
            markdown_result_file.write(my_format("From BAM:", gfa_inject_results["bam_inject"]))

        def extract_time(s):
            for l in s.decode('utf-8').splitlines():
                l = l.strip()
                if l.startswith("Elapsed"):  
                    timestamp = l.split(' ')[-1]
                   
                    # 0:00.11
                    t = 0.0
                    
                    if len(timestamp.split('.')) > 1:
                        t += float("0." + timestamp.split('.')[1])
                    timestamp = timestamp.split('.')[0]
                    # 0:00 
                    # Seconds per unit
                    base = 1
                    
                    for e in timestamp.split(':')[::-1]:
                        t += int(e)*base
                        base *= 60
                    return t

            print("No time line found.")  
            exit(1)

        def calc_total(obj):
            t = 0
            for _, v in obj.items():
                t += extract_time(v)
            return t
        
        markdown_result_file.write("# TOTAL WORKFLOW TIME\n\n\n")
        markdown_result_file.write(total_time_format("GBAM:", calc_total(gbam_results)))
        markdown_result_file.write(total_time_format("SAMBAMBA:", calc_total(sambamba_results)))
        markdown_result_file.write(total_time_format("SAMTOOLS:", calc_total(samtools_bam_results)))
        markdown_result_file.write(total_time_format("SAMTOOLS CRAM:", calc_total(samtools_cram_results)))
        
    print("Completed successfully.")
