"""
Author: Aaron Ho
Python Version: 3.7
"""
import argparse


def get_all_orfs(seq, minbp, nested):
    """
    Find and sort ORFs on both forward and reverse strands
    :param str seq: DNA sequence
    :param int minbp: Minimum bp length of ORF
    :param bool nested: True-Show Nested ORFs, False-Hide nested ORFs
    :return list: List of ORFs represented as tuples,
                    formatted as (START INDEX, END INDEX, FRAME, STRAND)
    """
    f_orfs = get_orfs(seq, minbp, nested)

    # Find reverse strand ORFs
    rc_seq = r_comp(seq)
    r_orfs_no_index = get_orfs(rc_seq, minbp, nested)
    r_orfs = [(len(rc_seq) + 1 - orf[1], len(rc_seq) + 1 - orf[0], orf[2], '-')
              for orf in r_orfs_no_index]   # Index ORFS in relation to forward strand
    r_orfs.sort(key=lambda orf: (orf[0], orf[1], orf[2], orf[3]))

    # Merge and sort both strands
    orf_list = [(orf[0], orf[1], orf[2], '+') for orf in f_orfs]
    orf_list.extend(r_orfs)
    orf_list.sort(key=lambda orf: (orf[0], orf[1], orf[2], orf[3]))

    return orf_list


def get_genes(seq, orf_list):
    """
    Retrieve gene sequences identified in ORF list
    :param str seq: Reference sequence used to make orf_list
    :param list orf_list: List of ORFs in
                            (START INDEX, END INDEX, FRAME, STRAND) format
    :return list: List of genes formatted as
                    (START INDEX, END INDEX, FRAME, STRAND, PROTEIN SEQUENCE)
    """
    gene_list = []

    for orf in orf_list:

        if orf[3] == '-':
            subseq = r_comp(seq[orf[0] - 1: orf[1]])
        else:
            subseq = seq[orf[0] - 1: orf[1]]

        protein_seq = nuc_to_pro(subseq)
        gene_list.append((orf[0], orf[1], orf[2], orf[3], protein_seq))

    return gene_list


def get_orfs(seq, minbp, nested):
    """
    :SOURCE: seq_tools/orf_finder.py
    """
    stop_c = ['TAG', 'TAA', 'TGA']
    orf_list = []

    for frame in range(3):
        start_pos, end_pos = ([] for l in range(2))

        can_start = True
        for i in range(frame, len(seq), 3):
            codon = seq[i:i + 3]

            if codon == 'ATG' and can_start:
                start_pos.append(i)
                can_start = False
            elif codon in stop_c:
                end_pos.append(i)

                can_start = True

        end_index = 0
        last_end = end_pos[len(end_pos) - 1]
        curr_orfs = []

        for i in start_pos:  # for every start position

            if i < end_pos[end_index]:  # if start less than end
                curr_orfs.append((i, end_pos[end_index]))
            elif i > last_end:
                break
            else:
                while end_pos[end_index] < i:
                    end_index += 1
                curr_orfs.append((i, end_pos[end_index]))

        prev_end = -1
        for i in curr_orfs:
            if (i[1] - i[0]) >= minbp:
                if nested:
                    orf_list.append((i[0] + 1, i[1], frame + 1))

                elif i[1] != prev_end:
                    orf_list.append((i[0] + 1, i[1], frame + 1))
            prev_end = i[1]

    orf_list.sort(key=lambda gene: (gene[0], gene[1], gene[2]))

    # This portion only runs if we want to remove nested
    if nested is False:
        orf_overlaps = sorted(orf_list,
                              key=lambda gene: (gene[1] - gene[0]), reverse=True)

        for i in orf_overlaps:
            orf_list[:] = [g for g in orf_list if g == i
                           or g[0] < i[0] or g[1] > i[1]]

    return orf_list


def nuc_to_pro(seq):
    """
    :SOURCE: seq_tools/dna_map.py
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


def r_comp(seq):
    """
    :SOURCE: seq_tools/dna_map.py
    """
    comp_base = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    seq_list = list(seq)
    seq_list.reverse()
    reverse = ''.join(seq_list)
    return ''.join([comp_base[i] for i in reverse])


def parse_fasta(fasta):
    """
    :SOURCE: file_tools/parse_fasta.py
    """
    with open(fasta) as file:
        seq = ''.join([line.strip()
                       for line in file if not line.startswith('>')])
    return seq


def main():
    parser = argparse.ArgumentParser(
        description='Identify genes on all strands and frames')
    parser.add_argument('-f', '--file', type=str, required=True,
                        help='Input DNA sequence file')
    parser.add_argument('-m', '--minbp', type=int, required=True,
                        help='Minimum ORF length')
    parser.add_argument('-n', '--nested', action='store_true',
                        help='Show nested ORFs')
    args = parser.parse_args()

    seq = parse_fasta(args.file)
    bp = args.minbp
    show_nested = args.nested

    orf_list = get_all_orfs(seq, bp, show_nested)
    gene_list = get_genes(seq, orf_list)

    print('Label\tStrand\tFrame\tStart\tStop\tSequence\n')

    g_index = 0
    for gene in gene_list:
        g_index += 1
        print(f'Gene {g_index}\t{gene[3]}\t{gene[2]}'
              f'\t{gene[0]}\t{gene[1]}\t{gene[4]}')


if __name__ == '__main__':
    main()
