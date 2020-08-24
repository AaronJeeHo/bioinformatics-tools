"""
Author: Aaron Ho
Python Version: 3.7
"""
import argparse


def parse_fastq(fastq):
    """
    Parse data from fastq
    :param str fastq: Fastq file
    :return dict: Dict containing Fastq records formatted as (SEQ, QUAL)
    """
    fq_records = {}

    with open(fastq) as file:
        l_count = 0
        curr_record = curr_seq = curr_qual = ''
        sep_found = False
        for line in file:
            l_count += 1
            if l_count % 4 == 1 and line.startswith('@'):
                curr_record = line[1:].rstrip()
            elif l_count % 4 == 2:
                curr_seq = line.rstrip()
            elif l_count % 4 == 3 and line.startswith('+') is True:
                sep_found = True
            elif l_count % 4 == 0 and sep_found is True:
                curr_qual = line.rstrip()
                fq_records[curr_record] = (curr_seq, curr_qual)
                curr_record = curr_seq = curr_qual = ''
                sep_found = False

    return fq_records


def main():
    parser = argparse.ArgumentParser(
        description='Parse sequences from FASTQ file')
    parser.add_argument('file', help='Input FASTQ sequence file')
    parser.add_argument('-f', '--fasta', action='store_true',
                        help='Output in FASTA format')
    args = parser.parse_args()

    records = parse_fastq(args.file)
    if args.fasta:
        for k in records.keys():
            print(f'>{k}\n{records[k][0]}\n')
    else:
        print('\n'.join(records[k][0] for k in records.keys()))


if __name__ == '__main__':
    main()
