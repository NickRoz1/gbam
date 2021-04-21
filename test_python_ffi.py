import gbam_tools

def convert(args):
    from gbam_tools import bam_to_gbam_python
    bam_to_gbam_python(str(args.input_file), str(args.output_file))
    
def read(args):
    from gbam_tools import Reader, ParsingTemplate

    from enum import IntEnum


    class Fields(IntEnum):
        POS = 1
        MAPQ = 3


    FIELDS_NUM = 17

    fields = [False] * FIELDS_NUM
    fields[Fields.POS] = True
    fields[Fields.MAPQ] = True

    parsing_template = ParsingTemplate(fields)

    reader = Reader("../test_data/res.gbam", parsing_template)
    pass

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

    subparsers = parser.add_subparsers(help='Choose one of the following options.', dest="command")
    
    parser_convert = subparsers.add_parser('conv', help='Convert BAM file to GBAM file.')
    parser_convert.add_argument("-i", "--input_file", help="BAM file to convert.", dest="input_file",
                                type=lambda x: is_valid_file(parser_convert, x), required=True)
    parser_convert.add_argument(
        "-o", "--output_file", type=str, help="GBAM file name to write out.", required=True)

    parser_read = subparsers.add_parser('read', help='Read GBAM file.')
    parser_read.add_argument("-i", "--input_file", help="GBAM file to read.", dest="input_file",
                             type=lambda x: is_valid_file(parser_convert, x), required=True)
    parser_read.add_argument('--no_mapq', dest='mapq', action='store_true',
                             help='Avoid reading MAPQ field.')
    parser_read.add_argument('--no_pos', dest='pos', action='store_true',
                             help='Avoid reading POS field.')
    import sys
    
    args = parser.parse_args(sys.argv[1:])
    # print(args)
    if args.command == 'conv':
        convert(args)
    elif args.command == 'read':
        args = parser_read.parse_args(sys.argv[1:])
        read(args)