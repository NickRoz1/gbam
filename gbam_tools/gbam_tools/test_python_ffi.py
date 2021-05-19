import gbam_tools
from tempfile import NamedTemporaryFile
import pysam
from enum import IntEnum
from array import array


# These values should have same values as in Rust version
class Fields(IntEnum):
    REFID = 0
    POS = 1
    MAPQ = 2 
    BIN = 3
    FLAGS = 4
    NEXTREFID = 5
    NEXTPOS = 6
    TLEN = 7
    READNAME = 8
    RAWCIGAR = 9
    RAWSEQUENCE = 10
    RAWQUAL = 11
    RAWTAGS = 12



FIELDS_NUM = 13

def convert(bam_path, gbam_path, compression = 'gzip'):
    from gbam_tools import bam_to_gbam_python
    bam_to_gbam_python(bam_path, gbam_path, compression)
    print("Conversion completed.")


def get_reader(path, parsing_tmplt):
    from gbam_tools import Reader, ParsingTemplate
    return Reader(path, parsing_tmplt)


# Converts BAM file to GBAM file, performs tests, deletes GBAM file.
def test(args):
    input_path = args.input_path

    output_file = NamedTemporaryFile()
    convert(input_path, output_file.name, 'lz4')
    output_file_name = output_file.name

    # test_combinations = 0

    # for field_to_omit in list(map(int, Fields)):
    # input_path = "../test_data/1.bam"
    # output_file_name = "../test_data/res.gbam"
    check_if_equal(input_path, output_file_name)

    print("Tests passed.")

def get_parsing_tmpl(fields_to_parse):
    from gbam_tools import ParsingTemplate
    fields = [False] * FIELDS_NUM
    for field in fields_to_parse:
        fields[field] = True
    return ParsingTemplate(fields)

# Checks if the data in both BAM and GBAM files is equal.


def check_if_equal(bam_path, gbam_path, no_check_fields=[]):
    # Suppress warnings to work with BAM files without index file. 
    # https://github.com/pysam-developers/pysam/issues/939#issuecomment-669016051
    save = pysam.set_verbosity(0)
    bam_file = pysam.AlignmentFile(bam_path, "rb")
    pysam.set_verbosity(save)

    fields_to_check = [field for field in list(map(int, Fields)) if field not in no_check_fields]

    gbam_file = get_reader(gbam_path, get_parsing_tmpl(fields_to_check)) 
    from gbam_tools import GbamRecord

    i = 0
    while True:
        cur_gbam = gbam_file.next_record()
        cur_bam = next(bam_file, None)
        if i > 0 and i % 100000 == 0:
            print('%d records are processed' % i)
        if cur_gbam == None or cur_bam == None:
            # Assert there is no records left
            assert(cur_gbam == cur_bam)
            break
        
        for field in fields_to_check:
            if field == Fields.REFID:
                assert(cur_bam.reference_id == cur_gbam.refid)
            if field == Fields.POS:
                assert(cur_bam.reference_start == cur_gbam.pos)
            if field == Fields.MAPQ:
                assert(cur_bam.mapping_quality == cur_gbam.mapq)
            if field == Fields.BIN:
                assert(cur_bam.bin == cur_gbam.bin)
            if field == Fields.FLAGS:
                assert(cur_bam.flag == cur_gbam.flag)
            if field == Fields.NEXTREFID:
                assert(cur_bam.next_reference_id == cur_gbam.next_ref_id)
            if field == Fields.NEXTPOS:
                assert(cur_bam.next_reference_start == cur_gbam.next_pos)
            if field == Fields.TLEN:
                assert(cur_bam.template_length == cur_gbam.tlen)
            if field == Fields.READNAME:
                assert(list(bytearray(cur_bam.query_name, 'utf8')) == cur_gbam.read_name[:-1])
            if field == Fields.RAWCIGAR:
                assert(cur_bam.cigarstring == cur_gbam.cigar)
            if field == Fields.RAWSEQUENCE:
                assert(cur_bam.query_sequence == cur_gbam.seq)
            if field == Fields.RAWQUAL:
                assert(cur_bam.query_qualities == array('B', cur_gbam.qual))
        i += 1
                
def is_valid_file(parser, arg):
    import os.path

    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)
    else:
        return arg


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Tool for creating and reading GBAM files",
    )

    subparsers = parser.add_subparsers(
        help='Choose one of the following options.', dest="command")
    subparsers.required = True

    parser_convert = subparsers.add_parser(
        'conv', help='Convert BAM file to GBAM file.')
    parser_convert.add_argument("-i", "--input_path", help="BAM file to convert.", dest="input_path",
                                type=lambda x: is_valid_file(parser_convert, x), required=True)
    parser_convert.add_argument(
        "-o", "--output_file", type=str, help="GBAM file name to write out.", required=True)
    parser_convert.add_argument(
        "--lz4", action='store_true', help="Use lz4 for compression")
    # parser_read = subparsers.add_parser('read', help='Read GBAM file.')
    # parser_read.add_argument("-i", "--input_path", help="GBAM file to read.", dest="input_path",
    #                          type=lambda x: is_valid_file(parser_convert, x), required=True)
    # parser_read.add_argument('--no_mapq', dest='mapq', action='store_true',
    #                          help='Avoid reading MAPQ field.')
    # parser_read.add_argument('--no_pos', dest='pos', action='store_true',
    #                          help='Avoid reading POS field.')

    parser_test = subparsers.add_parser('test', help='Test GBAM tools.')
    parser_test.add_argument("-i", "--input_path", help="GBAM file to perform tests on.", dest="input_path",
                             type=lambda x: is_valid_file(parser_convert, x), required=True)

    import sys

    args = parser.parse_args(sys.argv[1:])

    if args.command == 'conv':
        compression = 'gzip'
        if args.lz4:
            compression = 'lz4'
        convert(args.input_path, args.output_file, compression)
    # elif args.command == 'read':
    #     args = parser_read.parse_args(sys.argv[1:])
    #     read(args)
    elif args.command == 'test':
        test(args)
