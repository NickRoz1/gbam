# GBAM

This package includes tools for creating and manipulating GBAM files. GBAM is a new column oriented file format for storing binary alignment data. The way it stores information allows to skip unnecessary operations which are inevitable when using regular BAM files.

# Installation

## Install binary

```shell
git clone https://github.com/NickRoz1/gbam
cd gbam
cd gbam_binary
cargo build
cargo test
../target/debug/gbam_binary --help
```

You may need to set the C compiler with

```shell
env CC=gcc cargo build
```

## Build library

To build the `libgbam_tools.so` library you can run following command in the `gbam_tools` directory.

```bash
cargo build --release --features python-ffi
```



## With Python tools

```shell
git clone https://github.com/NickRoz1/gbam
cd gbam
cd gbam_tools
python3 -m venv env
source env/bin/activate
pip3 install maturin
maturin develop --release --cargo-extra-args="--features python-ffi"
```

You may need to install following packages:
```shell
sudo apt-get install liblzma-dev
sudo apt-get install libbz2-dev
```

# Usage

### To convert BAM file to GBAM file
```shell
python3 test_python_ffi.py conv -i test_data/input.bam -o test_data/output.gbam
```
### To test whether the GBAM file contains correct info and whether it iterates properly
```shell
python3 test_python_ffi.py test -i test_data/input.bam
```

# Development

## virtualenv

maturin can use python's virtualenv. See above instructions for Python install.

## GNU Guix

GNU Guix provides a full environment for development.  See
[.guix-build](./.guix-build)
