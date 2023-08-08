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


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Benchmark SAMTOOLS/SAMBAMBA/GBAM")
    parser.add_argument("--samtools_bin",  help="Path to samtools binary", required=False)
    parser.add_argument("--sambamba_bin",  help="Path to sambamba binary", required=False)
    parser.add_argument("--gbam_bin",  help="Path to gbam binary", required=True)

    parser.add_argument("--bam_file",  help="Path to BAM file to perform test on.", required=True)
    parser.add_argument("--result_dir",  help="Path to the directory where to save the resulting files", required=True)

    args = parser.parse_args()

    for p in [args.samtools_bin, args.sambamba_bin, args.gbam_bin, args.bam_file]:
        if not os.path.exists(args.result_dir):
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

    print("Completed successfully.")
        
    
