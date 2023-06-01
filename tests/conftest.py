from test_conversions import bam_file_path

def pytest_addoption(parser):
    parser.addoption("--bam-file-path", action="store", help="Provide BAM test file to perform test on. Currently only test_sort and test_bam_to_gbam_to_bam are using it. This option is needed to avoid running heavy fixture on big files.")
    parser.addoption(
        "--use-custom-file", action="store", default=bam_file_path, help="Provide a full path to a file you want to use to run all tests on."
    )