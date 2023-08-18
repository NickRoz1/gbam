# GBAM Benchmarks

Generated with:
```sh
tests/benchmark.py --gbam_bin target/release/gbam_binary --bam_file test_data/S288.3X.150.vs.pangenome.bam --result_dir benchmarking --samtools_bin /home/wrk/opt/samtools/bin/samtools --sambamba_bin /home/wrk/iwrk/opensource/code/D/sambamba/bin/sambamba-1.0.1-linux-amd64-static --gfainject_bin /home/wrk/iwrk/opensource/code/pangenome/gfainject/target/release/gfainject --gfa_file test_data/scerevisiae7.fa.gz.d1a145e.417fcdf.7493449.smooth.final.gfa --gbam_file test_data/S288.3X.150.vs.pangenome.gbam
```

# SORTING


## GBAM index-sort:
 ```
User time (seconds): 3.83
System time (seconds): 0.48
Percent of CPU this job got: 540%
Elapsed (wall clock) time (h:mm:ss or m:ss): 0:00.80
Maximum resident set size (kbytes): 485672
```

## SAMBAMBA sort:
 ```
User time (seconds): 13.35
System time (seconds): 0.54
Percent of CPU this job got: 752%
Elapsed (wall clock) time (h:mm:ss or m:ss): 0:01.84
Maximum resident set size (kbytes): 692012
```

## SAMTOOLS sort:
 ```
User time (seconds): 13.16
System time (seconds): 0.25
Percent of CPU this job got: 766%
Elapsed (wall clock) time (h:mm:ss or m:ss): 0:01.75
Maximum resident set size (kbytes): 669492
```

# FLAGSTAT


## GBAM flagstat:
 ```
User time (seconds): 0.08
System time (seconds): 0.02
Percent of CPU this job got: 216%
Elapsed (wall clock) time (h:mm:ss or m:ss): 0:00.04
Maximum resident set size (kbytes): 116536
```

## SAMBAMBA flagstat:
 ```
User time (seconds): 1.16
System time (seconds): 0.05
Percent of CPU this job got: 812%
Elapsed (wall clock) time (h:mm:ss or m:ss): 0:00.15
Maximum resident set size (kbytes): 11308
```

## SAMTOOLS flagstat:
 ```
User time (seconds): 1.77
System time (seconds): 0.11
Percent of CPU this job got: 817%
Elapsed (wall clock) time (h:mm:ss or m:ss): 0:00.23
Maximum resident set size (kbytes): 11564
```

# DEPTH


## GBAM depth:
 ```
User time (seconds): 0.94
System time (seconds): 0.64
Percent of CPU this job got: 109%
Elapsed (wall clock) time (h:mm:ss or m:ss): 0:01.44
Maximum resident set size (kbytes): 117968
```

## SAMBAMBA depth:
 ```
User time (seconds): 33.74
System time (seconds): 1.16
Percent of CPU this job got: 107%
Elapsed (wall clock) time (h:mm:ss or m:ss): 0:32.59
Maximum resident set size (kbytes): 24792
```

## SAMTOOLS depth:
 ```
User time (seconds): 4.42
System time (seconds): 0.71
Percent of CPU this job got: 183%
Elapsed (wall clock) time (h:mm:ss or m:ss): 0:02.79
Maximum resident set size (kbytes): 11420
```

# GFAINJECT


## From GBAM:
 ```
User time (seconds): 2.88
System time (seconds): 1.20
Percent of CPU this job got: 98%
Elapsed (wall clock) time (h:mm:ss or m:ss): 0:04.15
Maximum resident set size (kbytes): 125004
```

## From BAM:
 ```
User time (seconds): 4.97
System time (seconds): 1.31
Percent of CPU this job got: 99%
Elapsed (wall clock) time (h:mm:ss or m:ss): 0:06.29
Maximum resident set size (kbytes): 37512
```

