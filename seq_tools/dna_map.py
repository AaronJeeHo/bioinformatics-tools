"""
Author: Aaron Ho
Python Version: 3.7
"""
import argparse


def get_len(seq):
    """
    Measures sequence length
    :param str seq: DNA sequence
    :return int: Sequence length
    """
    return len(seq)


def nuc_count(seq):
    """
    Return nucleotide count in sequence
    :param str seq: DNA sequence
    :return list: List of nucleotide counts
    """
    seq_list = list(seq)
    return (seq_list.count('A'), seq_list.count('C'),
            seq_list.count('G'), seq_list.count('T'))


def to_rna(seq):
    """
    Convert DNA sequence to RNA
    :param str seq: DNA sequence
    :return str: RNA sequence
    """
    return seq.replace('T', 'U')


def r_comp(seq):
    """
    Returns reverse complementary sequence
    :param str seq: DNA sequence
    :return str: Reverse complementary sequence
    """
    comp_base = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    seq_list = list(seq)
    seq_list.reverse()
    reverse = ''.join(seq_list)
    return ''.join([comp_base[i] for i in reverse])


def nuc_to_pro(seq):
    """
    Converts DNA or RNA sequence to protein sequence
    :param str seq: DNA sequence w/ length divisible by 3
    :return str: Protein sequence
    """
    translate = {
        'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C',
        'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C',
        'TTA': 'L', 'TCA': 'S', 'TAA': '', 'TGA': '',
        'TTG': 'L', 'TCG': 'S', 'TAG': '', 'TGG': 'W',

        'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R',
        'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',
        'CTA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',
        'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',

        'ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT': 'S',
        'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',
        'ATA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',
        'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',

        'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G',
        'GTC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',
        'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',
        'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'
    }

    if len(seq) % 3 != 0:
        return ("Sequence not converted to codons, "
                "unable to group nucleotides in threes")

    if 'U' in seq:
        seq = seq.replace('U', 'T')

    return ''.join([translate[seq[i:i+3]] for i in range(0, len(seq), 3)])


def main():
    parser = argparse.ArgumentParser(description='Manipulate DNA Sequence')
    parser.add_argument('sequence', help='Input DNA sequence')

    # Options for sequence manipulation
    seq_action = parser.add_mutually_exclusive_group()
    seq_action.add_argument('-l', '--length', action='store_true',
                            help='Return sequence length')

    seq_action.add_argument('-n', '--nuc', action='store_true',
                            help='Return nucleotide counts')

    seq_action.add_argument('-r', '--rna', action='store_true',
                            help='Convert DNA to RNA')

    seq_action.add_argument('-c', '--comp', action='store_true',
                            help='Return reverse complementary strand')

    seq_action.add_argument('-p', '--protein', action='store_true',
                            help='Convert DNA to Protein')

    args = parser.parse_args()
    seq = args.sequence

    # Perform specified action
    if args.length:
        print(get_len(seq))
    elif args.rna:
        print(to_rna(seq))
    elif args.comp:
        print(r_comp(seq))
    elif args.protein:
        print(nuc_to_pro(seq))
    elif args.nuc:
        nuc_list = nuc_count(seq)
        print(f'A: {nuc_list[0]}')
        print(f'C: {nuc_list[1]}')
        print(f'G: {nuc_list[2]}')
        print(f'T: {nuc_list[3]}')
    else:
        print("Please Input Command")


if __name__ == '__main__':
    main()
