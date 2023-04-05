from pathlib import Path
import subprocess
from tempfile import NamedTemporaryFile

cur_file_path = Path(__file__).parent.absolute()

test_data_folder = cur_file_path.parent/"test_data"
binary_path = cur_file_path.parent/"target"/"release"/"gbam_binary"

bam_file_path = test_data_folder/"testing_gbam_to_bam"/"wgEncodeUwRepliSeqGm12878G1bAlnRep1.bam"

def test_bam_to_gbam_to_bam():
    gbam_file = NamedTemporaryFile()
    bam_file_from_gbam = NamedTemporaryFile(suffix=".bam")

    subprocess.run([binary_path, bam_file_path, "-c", gbam_file.name]) 
    subprocess.run([binary_path, "--convert-to-bam", gbam_file.name, bam_file_from_gbam.name]) 

    view_of_original = subprocess.check_output(["samtools", "view", str(bam_file_path)], stderr=subprocess.STDOUT)
    view_of_result = subprocess.check_output(["samtools", "view", bam_file_from_gbam.name], stderr=subprocess.STDOUT)

    assert(len(view_of_original) > 0)
    assert(view_of_original == view_of_result)

def test_flagstat():
    flagstat_bam = test_data_folder/"wgEncodeRikenCageHuvecCellPapAlnRep1.bam"

    gbam_file = NamedTemporaryFile()

    subprocess.run([binary_path, flagstat_bam, "-c", gbam_file.name]) 
    print(gbam_file.name)

    view_of_original = subprocess.check_output(["samtools", "flagstat", str(flagstat_bam)], stderr=subprocess.STDOUT)
    view_of_result = subprocess.check_output([binary_path, "--flagstat", gbam_file.name], stderr=subprocess.STDOUT)

    assert(len(view_of_original) > 0)
    assert(view_of_original == view_of_result)