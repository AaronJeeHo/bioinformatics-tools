"""
Author: Aaron Ho
Python Version: 3.7
"""
import argparse


def parse_fasta(fasta):
    """
    Parse Fasta File
    :param str fasta: File Name
    :return str: Sequence
    """
    with open(fasta) as file:
        seq = ''.join([line.strip()
                       for line in file if not line.startswith('>')])
    return seq


def parse_multiseq(fasta):
    """
    Parse Fasta containing multiple sequences denoted by '>'
    :param str fasta: File Name
    :return list: List containing sequences
    """
    seq_list = []
    curr_seq = None
    with open(fasta) as file:
        for line in file:
            if line.startswith('>'):
                if curr_seq is not None:
                    seq_list.append(''.join(curr_seq))
                curr_seq = []
            else:
                curr_seq.append(line.strip())

    if curr_seq is not None:
        seq_list.append(''.join(curr_seq))
    return seq_list


def main():
    parser = argparse.ArgumentParser(
        description='Extract sequences from FASTA file')
    parser.add_argument('file', help='Input FASTA sequence file')
    parser.add_argument('-m', '--multi', action='store_true',
                        help='File contains multiple sequences')
    args = parser.parse_args()

    if args.multi:
        print('\n'.join(parse_multiseq(args.file)))
    else:
        print(parse_fasta(args.file))


if __name__ == '__main__':
    main()
