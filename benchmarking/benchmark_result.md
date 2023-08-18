# GBAM Benchmarks

Generated with:
```sh
tests/benchmark.py --gbam_bin target/release/gbam_binary --bam_file test_data/S288.3X.150.vs.pangenome.bam --result_dir benchmarking --samtools_bin /home/wrk/opt/samtools/bin/samtools --sambamba_bin /home/wrk/iwrk/opensource/code/D/sambamba/bin/sambamba-1.0.1-linux-amd64-static --gfainject_bin /home/wrk/iwrk/opensource/code/pangenome/gfainject/target/release/gfainject --gfa_file test_data/scerevisiae7.fa.gz.d1a145e.417fcdf.7493449.smooth.final.gfa --gbam_file test_data/S288.3X.150.vs.pangenome.gbam
```

# SORTING


## GBAM index-sort:
 ```
User time (seconds): 3.75
System time (seconds): 0.49
Percent of CPU this job got: 553%
Elapsed (wall clock) time (h:mm:ss or m:ss): 0:00.76
Maximum resident set size (kbytes): 485596
```

## SAMBAMBA sort:
 ```
User time (seconds): 13.50
System time (seconds): 0.49
Percent of CPU this job got: 766%
Elapsed (wall clock) time (h:mm:ss or m:ss): 0:01.82
Maximum resident set size (kbytes): 686640
```

## SAMTOOLS sort:
 ```
User time (seconds): 13.50
System time (seconds): 0.26
Percent of CPU this job got: 767%
Elapsed (wall clock) time (h:mm:ss or m:ss): 0:01.79
Maximum resident set size (kbytes): 670796
```

# FLAGSTAT


## GBAM flagstat:
 ```
User time (seconds): 0.08
System time (seconds): 0.02
Percent of CPU this job got: 211%
Elapsed (wall clock) time (h:mm:ss or m:ss): 0:00.05
Maximum resident set size (kbytes): 115576
```

## SAMBAMBA flagstat:
 ```
User time (seconds): 1.23
System time (seconds): 0.04
Percent of CPU this job got: 812%
Elapsed (wall clock) time (h:mm:ss or m:ss): 0:00.15
Maximum resident set size (kbytes): 11308
```

## SAMTOOLS flagstat:
 ```
User time (seconds): 1.60
System time (seconds): 0.09
Percent of CPU this job got: 803%
Elapsed (wall clock) time (h:mm:ss or m:ss): 0:00.21
Maximum resident set size (kbytes): 11496
```

# DEPTH


## GBAM depth:
 ```
User time (seconds): 0.93
System time (seconds): 0.64
Percent of CPU this job got: 110%
Elapsed (wall clock) time (h:mm:ss or m:ss): 0:01.43
Maximum resident set size (kbytes): 115380
```

## SAMBAMBA depth:
 ```
User time (seconds): 33.74
System time (seconds): 1.12
Percent of CPU this job got: 107%
Elapsed (wall clock) time (h:mm:ss or m:ss): 0:32.54
Maximum resident set size (kbytes): 26884
```

## SAMTOOLS depth:
 ```
User time (seconds): 4.40
System time (seconds): 0.71
Percent of CPU this job got: 179%
Elapsed (wall clock) time (h:mm:ss or m:ss): 0:02.85
Maximum resident set size (kbytes): 11460
```

# GFAINJECT


## From GBAM:
 ```
User time (seconds): 5.02
System time (seconds): 1.27
Percent of CPU this job got: 99%
Elapsed (wall clock) time (h:mm:ss or m:ss): 0:06.30
Maximum resident set size (kbytes): 37428
```

## From BAM:
 ```
User time (seconds): 2.90
System time (seconds): 1.11
Percent of CPU this job got: 98%
Elapsed (wall clock) time (h:mm:ss or m:ss): 0:04.07
Maximum resident set size (kbytes): 125940
```

