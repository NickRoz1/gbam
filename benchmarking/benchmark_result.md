# GBAM Benchmarks

Generated with:
```sh
tests/benchmark.py --gbam_bin target/release/gbam_binary --bam_file test_data/S288.3X.150.vs.pangenome.bam --result_dir benchmarking --samtools_bin /home/wrk/opt/samtools/bin/samtools --sambamba_bin /home/wrk/iwrk/opensource/code/D/sambamba/bin/sambamba-1.0.1-linux-amd64-static --gfainject_bin /home/wrk/iwrk/opensource/code/pangenome/gfainject/target/release/gfainject --gfa_file test_data/scerevisiae7.fa.gz.d1a145e.417fcdf.7493449.smooth.final.gfa --gbam_file test_data/S288.3X.150.vs.pangenome.gbam
```

# SORTING


## GBAM index-sort:
 ```
User time (seconds): 3.69
System time (seconds): 0.44
Percent of CPU this job got: 519%
Elapsed (wall clock) time (h:mm:ss or m:ss): 0:00.79
Maximum resident set size (kbytes): 488104
```

## SAMBAMBA sort:
 ```
User time (seconds): 13.31
System time (seconds): 0.44
Percent of CPU this job got: 752%
Elapsed (wall clock) time (h:mm:ss or m:ss): 0:01.82
Maximum resident set size (kbytes): 685068
```

## SAMTOOLS sort:
 ```
User time (seconds): 13.74
System time (seconds): 0.30
Percent of CPU this job got: 776%
Elapsed (wall clock) time (h:mm:ss or m:ss): 0:01.80
Maximum resident set size (kbytes): 667852
```

# FLAGSTAT


## GBAM flagstat:
 ```
User time (seconds): 0.07
System time (seconds): 0.02
Percent of CPU this job got: 267%
Elapsed (wall clock) time (h:mm:ss or m:ss): 0:00.03
Maximum resident set size (kbytes): 99676
```

## SAMBAMBA flagstat:
 ```
User time (seconds): 1.29
System time (seconds): 0.03
Percent of CPU this job got: 852%
Elapsed (wall clock) time (h:mm:ss or m:ss): 0:00.15
Maximum resident set size (kbytes): 11224
```

## SAMTOOLS flagstat:
 ```
User time (seconds): 1.55
System time (seconds): 0.06
Percent of CPU this job got: 880%
Elapsed (wall clock) time (h:mm:ss or m:ss): 0:00.18
Maximum resident set size (kbytes): 11608
```

# DEPTH


## GBAM depth:
 ```
User time (seconds): 1.02
System time (seconds): 0.53
Percent of CPU this job got: 112%
Elapsed (wall clock) time (h:mm:ss or m:ss): 0:01.38
Maximum resident set size (kbytes): 115264
```

## SAMBAMBA depth:
 ```
User time (seconds): 33.90
System time (seconds): 1.05
Percent of CPU this job got: 107%
Elapsed (wall clock) time (h:mm:ss or m:ss): 0:32.59
Maximum resident set size (kbytes): 27388
```

## SAMTOOLS depth:
 ```
User time (seconds): 4.51
System time (seconds): 0.68
Percent of CPU this job got: 173%
Elapsed (wall clock) time (h:mm:ss or m:ss): 0:02.99
Maximum resident set size (kbytes): 11448
```

# GFAINJECT


## From GBAM:
 ```
User time (seconds): 4.95
System time (seconds): 1.36
Percent of CPU this job got: 99%
Elapsed (wall clock) time (h:mm:ss or m:ss): 0:06.32
Maximum resident set size (kbytes): 39456
```

## From BAM:
 ```
User time (seconds): 2.87
System time (seconds): 1.17
Percent of CPU this job got: 99%
Elapsed (wall clock) time (h:mm:ss or m:ss): 0:04.05
Maximum resident set size (kbytes): 123524
```

