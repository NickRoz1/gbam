# GBAM

This package includes tools for creating and manipulating GBAM files. GBAM is a new column oriented file format for storing binary alignment data. The way it stores information allows to skip unnecessary operations which are inevitable when using regular BAM files.

# Installation

## Install binary

```shell
git clone https://github.com/NickRoz1/gbam
cd gbam
cargo build --release
./target/release/gbam_binary --help
```

You may need to install following packages:
```shell
sudo apt-get install liblzma-dev
sudo apt-get install libbz2-dev
```

# Usage

### To convert BAM file to GBAM file
```shell
# Simply convert
time ./target/release/gbam_binary -c test_data/wgEncodeUwRepliSeqGm12878G1bAlnRep1.bam -o test_data/res.gbam

#Sort before writing (sort by Name and match mates)
time ./target/release/gbam_binary -c -s test_data/wgEncodeUwRepliSeqGm12878G1bAlnRep1.bam -o test_data/res.gbam
```
### To run pytests, just run this command
```shell
pytest
```

# Development

## GNU Guix - this is not tested.

GNU Guix provides a full environment for development.  See
[.guix-build](./.guix-build)
