import argparse
import subprocess
import os
from pathlib import Path
import time

samtools_results = {}
sambamba_results = {}
gbam_results = {}
gfa_inject_results = {}

keyword = "Elapsed"
my_format = lambda name, out: f"## {name}\n ```\n{out.decode('utf-8')}```\n\n"
my_format_only_elapsed = lambda name, out: f"## {name}\n ```\n{[s.strip() for s in out.decode('utf-8').splitlines() if keyword in s][0]}\n```\n\n"

def gbam(bin_path, bam_path, res_path):
    start = time.time()

    gbam_path = res_path/"index_sort_out.gbam"
    depth_gbam_out = res_path/"gbam_depth_out.txt"
    
    gbam_results["index-sort"]  = subprocess.check_output([f"/usr/bin/time -v {bin_path} {bam_path} -c -s --sort-temp-mode ram -o {gbam_path} --index-sort"], shell=True, stderr=subprocess.STDOUT)
    gbam_results["flagstat"]  = subprocess.check_output([f"/usr/bin/time -v {bin_path} --flagstat {gbam_path}"], shell=True, stderr=subprocess.STDOUT)
    gbam_results["depth"] = subprocess.check_output([f"/usr/bin/time -v {bin_path} --depth {gbam_path} --thread-num 12 --index-file {gbam_path.with_suffix('.gbam.gbai')} > {depth_gbam_out}"], shell=True, stderr=subprocess.STDOUT)
    
    print(f"Completed GBAM benchmarking, took: {time.time()-start} seconds")

def samtools(bin_path, bam_path, res_path):
    start = time.time()

    samtools_path = res_path/"sort_out.samtools.bam"
    depth_samtools_out = res_path/"samtools_depth_out.txt"
    
    samtools_results["sort"] = subprocess.check_output([f"/usr/bin/time -v {bin_path} sort -@ 8 {bam_path} -o {samtools_path}"], shell=True, stderr=subprocess.STDOUT)
    samtools_results["flagstat"] = subprocess.check_output([f"/usr/bin/time -v {bin_path} flagstat -@ 8 {samtools_path}"], shell=True, stderr=subprocess.STDOUT)
    samtools_results["depth"] = subprocess.check_output([f"/usr/bin/time -v {bin_path} depth -@ 8 {samtools_path} > {depth_samtools_out}"], shell=True, stderr=subprocess.STDOUT)
    
    print(f"Completed SAMTOOLS benchmarking, took: {time.time()-start} seconds")

def sambamba(bin_path, bam_path, res_path):
    start = time.time()

    sambamba_path = res_path/"sort_out.sambamba.bam"
    depth_sambamba_out = res_path/"sambamba_depth_out.txt"
    
    sambamba_results["sort"] = subprocess.check_output([f"/usr/bin/time -v {bin_path} sort -t 8 -m 90GB {bam_path} -o {sambamba_path}"], shell=True, stderr=subprocess.STDOUT)
    sambamba_results["flagstat"] = subprocess.check_output([f"/usr/bin/time -v {bin_path} flagstat -t 8 {sambamba_path}"], shell=True, stderr=subprocess.STDOUT)
    sambamba_results["depth"] = subprocess.check_output([f"/usr/bin/time -v {bin_path} depth base -t 8 {sambamba_path} > {depth_sambamba_out}"], shell=True, stderr=subprocess.STDOUT)
    
    print(f"Completed SAMBAMBA benchmarking, took: {time.time()-start} seconds")

def gfainject(bin_path, gbam_path, gfa_path, bam_path, res_path):
    gfainject_from_bam_path = res_path/"injection_from_bam.gfainject.gaf"
    start = time.time()
    gfa_inject_results["bam_inject"] = subprocess.check_output([f"/usr/bin/time -v {bin_path} --gfa {gfa_path} --bam {bam_path }> {gfainject_from_bam_path}"], shell=True, stderr=subprocess.STDOUT)
    print(f"Completed injection from BAM, took: {time.time()-start} seconds")

    gfainject_from_gbam_path = res_path/"injection_from_gbam.gfainject.gaf"
    start = time.time()
    gfa_inject_results["gbam_inject"] = subprocess.check_output([f"/usr/bin/time -v {bin_path} --gfa {gfa_path} --gbam {gbam_path }> {gfainject_from_gbam_path}"], shell=True, stderr=subprocess.STDOUT)
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

    args = parser.parse_args()

    for p in [args.samtools_bin, args.sambamba_bin, args.gbam_bin, args.bam_file, args.gfainject_bin, args.gbam_file, args.gfa_file]: 
        print(p)
        if not os.path.exists(p):
            print(f"Path {p} does not exist, terminating.")
            exit()

    if not os.path.exists(args.result_dir):
        os.makedirs(args.result_dir)
        print(f"Directory '{args.result_dir}' created.")

    with open(Path(args.result_dir)/"benchmark_result.md", 'w') as markdown_result_file:
        same_params = (args.bam_file, Path(args.result_dir))
        gbam(args.gbam_bin, *same_params)
        samtools(args.samtools_bin, *same_params)
        sambamba(args.sambamba_bin, *same_params)
        gfainject(args.gfainject_bin, args.gbam_file, args.gfa_file, *same_params)

        # Print resulting file
        markdown_result_file.write("# SORTING\n\n\n")
        markdown_result_file.write(my_format_only_elapsed("GBAM index-sort:", gbam_results["index-sort"]))
        markdown_result_file.write(my_format_only_elapsed("SAMBAMBA sort:", sambamba_results["sort"]))
        markdown_result_file.write(my_format_only_elapsed("SAMTOOLS sort:", samtools_results["sort"]))

        markdown_result_file.write("# FLAGSTAT\n\n\n")
        markdown_result_file.write(my_format_only_elapsed("GBAM flagstat:", gbam_results["flagstat"]))
        markdown_result_file.write(my_format_only_elapsed("SAMBAMBA flagstat:", sambamba_results["flagstat"]))
        markdown_result_file.write(my_format_only_elapsed("SAMTOOLS flagstat:", samtools_results["flagstat"]))

        markdown_result_file.write("# DEPTH\n\n\n")
        markdown_result_file.write(my_format_only_elapsed("GBAM depth:", gbam_results["depth"]))
        markdown_result_file.write(my_format_only_elapsed("SAMBAMBA depth:", sambamba_results["depth"]))
        markdown_result_file.write(my_format_only_elapsed("SAMTOOLS depth:", samtools_results["depth"]))

        markdown_result_file.write("# GFAINJECT\n\n\n")
        markdown_result_file.write(my_format_only_elapsed("From GBAM:", gfa_inject_results["bam_inject"]))
        markdown_result_file.write(my_format_only_elapsed("From BAM:", gfa_inject_results["gbam_inject"]))

    print("Completed successfully.")
        
    
