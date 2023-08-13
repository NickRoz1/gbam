import argparse
import subprocess
import os
from pathlib import Path
import time

def gbam(bin_path, bam_path, res_path, md_file):
    start = time.time()

    gbam_path = res_path/"index_sort_out.gbam"
    depth_gbam_out = res_path/"gbam_depth_out.txt"
    
    md_file.write("# GBAM\n\n\n")
    my_format = lambda name, out: f"## {name}\n ```\n{out.decode('utf-8')}```\n\n"

    output = subprocess.check_output([f"/usr/bin/time -v {bin_path} {bam_path} -c -s --sort-temp-mode ram -o {gbam_path} --index-sort"], shell=True, stderr=subprocess.STDOUT)
    md_file.write(my_format("Converting BAM to GBAM via index-sort:", output))
    output = subprocess.check_output([f"/usr/bin/time -v {bin_path} --flagstat {gbam_path}"], shell=True, stderr=subprocess.STDOUT)
    md_file.write(my_format("GBAM flagstat:", output))
    output = subprocess.check_output([f"/usr/bin/time -v {bin_path} --depth {gbam_path} --thread-num 12 --index-file {gbam_path.with_suffix('.gbam.gbai')} > {depth_gbam_out}"], shell=True, stderr=subprocess.STDOUT)
    md_file.write(my_format("GBAM depth:", output))
    
    print(f"Completed GBAM benchmarking, took: {time.time()-start} seconds")

def samtools(bin_path, bam_path, res_path, md_file):
    start = time.time()

    samtools_path = res_path/"sort_out.samtools.bam"
    depth_samtools_out = res_path/"samtools_depth_out.txt"
    
    md_file.write("# SAMTOOLS\n\n\n")
    my_format = lambda name, out: f"## {name}\n ```\n{out.decode('utf-8')}```\n\n"

    output = subprocess.check_output([f"/usr/bin/time -v {bin_path} sort -@ 8 {bam_path} -o {samtools_path}"], shell=True, stderr=subprocess.STDOUT)
    md_file.write(my_format("Sorting BAM:", output))
    output = subprocess.check_output([f"/usr/bin/time -v {bin_path} flagstat -@ 8 {samtools_path}"], shell=True, stderr=subprocess.STDOUT)
    md_file.write(my_format("Flagstat:", output))
    output = subprocess.check_output([f"/usr/bin/time -v {bin_path} depth -@ 8 {samtools_path} > {depth_samtools_out}"], shell=True, stderr=subprocess.STDOUT)
    md_file.write(my_format("Depth:", output))
    
    print(f"Completed SAMTOOLS benchmarking, took: {time.time()-start} seconds")

def sambamba(bin_path, bam_path, res_path, md_file):
    start = time.time()

    sambamba_path = res_path/"sort_out.sambamba.bam"
    depth_sambamba_out = res_path/"sambamba_depth_out.txt"
    
    md_file.write("# SAMBAMBA\n\n\n")
    my_format = lambda name, out: f"## {name}\n ```\n{out.decode('utf-8')}```\n\n"

    output = subprocess.check_output([f"/usr/bin/time -v {bin_path} sort -t 8 -m 90GB {bam_path} -o {sambamba_path}"], shell=True, stderr=subprocess.STDOUT)
    md_file.write(my_format("Sorting BAM:", output))
    output = subprocess.check_output([f"/usr/bin/time -v {bin_path} flagstat -t 8 {sambamba_path}"], shell=True, stderr=subprocess.STDOUT)
    md_file.write(my_format("Flagstat:", output))
    output = subprocess.check_output([f"/usr/bin/time -v {bin_path} depth base -t 8 {sambamba_path} > {depth_sambamba_out}"], shell=True, stderr=subprocess.STDOUT)
    md_file.write(my_format("Depth:", output))
    
    print(f"Completed SAMTOOLS benchmarking, took: {time.time()-start} seconds")

def gfainject(bin_path, gbam_path, gfa_path, bam_path, res_path, md_file):
    my_format = lambda name, out: f"## {name}\n ```\n{out.decode('utf-8')}```\n\n"
    
    gfainject_from_bam_path = res_path/"injection_from_bam.gfainject.gaf"
    md_file.write("# GFAINJECT from BAM\n\n\n")
    start = time.time()
    output = subprocess.check_output([f"/usr/bin/time -v {bin_path} --gfa {gfa_path} --bam {bam_path }> {gfainject_from_bam_path}"], shell=True, stderr=subprocess.STDOUT)
    print(f"Completed injection from BAM, took: {time.time()-start} seconds")
    md_file.write(my_format("Injection from BAM:", output))   

    gfainject_from_gbam_path = res_path/"injection_from_gbam.gfainject.gaf"
    md_file.write("# GFAINJECT from GBAM\n\n\n")
    start = time.time()
    output = subprocess.check_output([f"/usr/bin/time -v {bin_path} --gfa {gfa_path} --gbam {gbam_path }> {gfainject_from_gbam_path}"], shell=True, stderr=subprocess.STDOUT)
    print(f"Completed injection from GBAM, took: {time.time()-start} seconds")
    md_file.write(my_format("Injection from GBAM:", output))   

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

    for p in [args.samtools_bin, args.sambamba_bin, args.gfainject_bin, args.gbam_bin, args.bam_file, args.gbam_file, args.gfa_file]:
        if not os.path.exists(p):
            print(f"Path {p} does not exist, terminating.")
            exit()

    if not os.path.exists(args.result_dir):
        os.makedirs(args.result_dir)
        print(f"Directory '{args.result_dir}' created.")

    with open(Path(args.result_dir)/"benchmark_result.md", 'w') as markdown_result_file:
        same_params = (args.bam_file, Path(args.result_dir), markdown_result_file)
        gbam(args.gbam_bin, *same_params)
        samtools(args.samtools_bin, *same_params)
        sambamba(args.sambamba_bin, *same_params)
        gfainject(args.gfainject_bin, args.gbam_file, args.gfa_file, *same_params)

    print("Completed successfully.")
        
    
