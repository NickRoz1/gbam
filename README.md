# GBAM

This package includes tools for creating and manipulating GBAM files. GBAM is a new column oriented file format for storing binary alignment data. The way it stores information allows to skip unnecessary operations which are inevitable when using regular BAM files.

# Installation

```shell
git clone https://github.com/NickRoz1/gbam
cd gbam
cd gbam_tools
python3 -m venv env
source env/bin/activate
pip3 install maturin
maturin develop --release
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