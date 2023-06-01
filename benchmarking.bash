#!/usr/bin/bash
#### Provide a path to a file.
#### To download a file use: 
# wget URL

input_path=$1
output_path=$2
gbam_binary_path=$3
sambamba_binary_path=$4

# Perform operations using input and output paths
echo "Input path: $input_path"
echo "Output path (should be without format suffix): $output_path"
echo "GBAM binary path: $gbam_binary_path"
echo "Sambamba binary path: $sambamba_binary_path"
echo "Temp dir is set as: $TMPDIR"

samtools_sorted_bam_file="$output_path.sorted.samtools.bam"
sambamba_sorted_bam_file="$output_path.sorted.sambamba.bam"
gbam_file="$output_path.gbam"
gbam_index_file="$output_path.gbam.gbai"

samtools_depth_file="$output_path.samtools.depth"
sambamba_depth_file="$output_path.sambamba.depth"
gbam_depth_file="$output_path.gbam.depth"

echo "---------------------------------------------------------------------------------------------------------------"
echo ""
echo "Benchmark gbam index-sort versus samtools and sambamba sort (you have to be built on feature/index_sort branch)..."
echo ""
echo ""
echo "---------------------------------------------------------------------------------------------------------------"

/usr/bin/time -v $gbam_binary_path $input_path -c -s --sort-temp-mode ram -o $gbam_file --index-sort
echo "---------------------------------------------------------------------------------------------------------------"
/usr/bin/time -v samtools sort -@ 8 $input_path -o $samtools_sorted_bam_file
echo "---------------------------------------------------------------------------------------------------------------"
/usr/bin/time -v $sambamba_binary_path sort -t 8 $input_path -o $sambamba_sorted_bam_file -m 90GB

echo "---------------------------------------------------------------------------------------------------------------"
echo ""
echo "Benchmark gbam flagstat versus samtools and sambamba flagstat..."
echo ""
echo ""
echo "---------------------------------------------------------------------------------------------------------------"

/usr/bin/time -v $gbam_binary_path --flagstat $gbam_file
echo "---------------------------------------------------------------------------------------------------------------"
/usr/bin/time -v samtools flagstat -@ 8 $samtools_sorted_bam_file
echo "---------------------------------------------------------------------------------------------------------------"
/usr/bin/time -v $sambamba_binary_path flagstat -t 8 $sambamba_sorted_bam_file

echo "---------------------------------------------------------------------------------------------------------------"
echo ""
echo "Benchmark gbam depth versus samtools and sambamba depth... (sambamba is disabled by default due to taking too long to complete)"
echo ""
echo ""
echo "---------------------------------------------------------------------------------------------------------------"

/usr/bin/time -v $gbam_binary_path --depth $gbam_file --thread-num 12 --index-file $gbam_index_file > $gbam_depth_file
echo "---------------------------------------------------------------------------------------------------------------"
/usr/bin/time -v samtools depth -@ 8 $samtools_sorted_bam_file > $samtools_depth_file
echo "---------------------------------------------------------------------------------------------------------------"
/usr/bin/time -v $sambamba_binary_path depth base -t 8 $sambamba_sorted_bam_file > $sambamba_depth_file
