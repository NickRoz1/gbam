from pathlib import Path
import subprocess
from tempfile import NamedTemporaryFile, TemporaryDirectory
import pytest 
import gzip

cur_file_path = Path(__file__).parent.absolute()

test_data_folder = cur_file_path.parent/"test_data"
binary_path = cur_file_path.parent/"target"/"release"/"gbam_binary"

bam_file_path = test_data_folder/"testing_gbam_to_bam"/"wgEncodeUwRepliSeqGm12878G1bAlnRep1.bam"
bam_file_sorted_path = test_data_folder/"wgEncodeRikenCageHuvecCellPapAlnRep1.sorted.bam"

gbam_file = None
gbam_file_sorted = None
@pytest.fixture(scope='module', autouse=True)
def gen_gbam_file():
    global gbam_file, gbam_file_sorted
    gbam_file = NamedTemporaryFile()
    subprocess.run([binary_path, bam_file_path, "-c", "-o", gbam_file.name]) 
    gbam_file_sorted = NamedTemporaryFile()
    subprocess.run([binary_path, bam_file_sorted_path, "-c", "-o", gbam_file_sorted.name]) 

def test_bam_to_gbam_to_bam():
    bam_file_from_gbam = NamedTemporaryFile(suffix=".bam")

    subprocess.run([binary_path, "--convert-to-bam", gbam_file.name, "-o", bam_file_from_gbam.name]) 

    view_of_original = subprocess.check_output(["samtools", "view", str(bam_file_path)], stderr=subprocess.STDOUT)
    view_of_result = subprocess.check_output(["samtools", "view", bam_file_from_gbam.name], stderr=subprocess.STDOUT)

    assert(len(view_of_original) > 0)
    assert(view_of_original == view_of_result)

def test_flagstat():
    view_of_original = subprocess.check_output(["samtools", "flagstat", str(bam_file_path)], stderr=subprocess.STDOUT)
    view_of_result = subprocess.check_output([binary_path, "--flagstat", gbam_file.name], stderr=subprocess.STDOUT)

    assert(len(view_of_original) > 0)
    assert(view_of_original == view_of_result)


# Testing against mosdepth.
def test_depth():    
    temp_dir = TemporaryDirectory()
    dir_path = Path(temp_dir.name)
    mosdepth_prefix = "testing_depth_pytest"
    mosdepth_suffix = ".per-base.bed.gz"
    mosdepth_file = (dir_path/(mosdepth_prefix+mosdepth_suffix)).as_posix()

    subprocess.check_call(['mosdepth', '-x', mosdepth_prefix, bam_file_sorted_path], cwd=temp_dir.name)
    bed_gz = NamedTemporaryFile()
    subprocess.check_call([binary_path, gbam_file_sorted.name, '-d', '-o', bed_gz.name])
    
    with gzip.open(bed_gz, "rb") as gbam_res, gzip.open(mosdepth_file, "rb") as mosdepth_res:
        # Read the decompressed data from both files
        decompressed_gbam_res = gbam_res.read()
        decompressed_mosdepth_res= mosdepth_res.read()

    assert(len(decompressed_gbam_res) > 0)
    assert(decompressed_gbam_res == decompressed_mosdepth_res)