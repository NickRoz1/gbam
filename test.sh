source gbam_tools/env/bin/activate
cd gbam_tools
maturin develop --release --cargo-extra-args="--features python-ffi"
cd gbam_tools
python3 test_python_ffi.py ../../test_data/wgEncodeUwRepliSeqGm12878G1bAlnRep1.bam ../../test_data/wgEncodeUwRepliSeqGm12878G1bAlnRep1.sorted.bam