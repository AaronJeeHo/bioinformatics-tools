"""
Author: Aaron Ho
Python Version: 3.7
"""
import argparse


class AlignmentGraph:
    """
    Class that creates dynamic programming graph
    and handles alignment scoring through
    Needleman-Wunsch algorithm
    """
    def __init__(self, col, row, match, mis, indel):
        """
        Create and score alignment graph
        :param string col: DNA Sequence
        :param string row: DNA Sequence
        :param int match: Score value for match
        :param int mis: Score penalty for mismatch
        :param int indel: Score penalty for indel
        """
        self.matrix = []
        self.col = col
        self.row = row
        self.match = match
        self.mis = mis
        self.indel = indel
        self.source = (0, 0)
        self.sink = (len(col), len(row))

        # fill matrix with empty rows
        for j in range(len(row) + 1):
            self.matrix.append([])

        # start scoring
        for j in range(len(row) + 1):
            curr_row = self.matrix[j]

            if j == 0:
                curr_row.append(0)
                for i in range(1, len(col) + 1):
                    curr_row.append(self.get_node(i - 1, j) + indel)

            else:
                for i in range(len(col) + 1):
                    if i == 0:
                        curr_row.append(self.get_node(i, j - 1) + indel)
                    else:
                        scores = []  # [diag , left , up]

                        d_score = self.get_node(i - 1, j - 1)
                        if col[i - 1] == row[j - 1]:
                            scores.append(d_score + match)
                        else:
                            scores.append(d_score + mis)

                        scores.append(self.get_node(i, j - 1) + indel)
                        scores.append(self.get_node(i - 1, j) + indel)

                        curr_row.append(max(scores))

    def back_trace(self):
        """
        Backtrack and find alignment path
        :return list: List containing [SCORE, ALIGNMENT LEN,
                        COLUMN ALIGNMENT, ROW ALIGNMENT]
        """
        col_seq = []
        row_seq = []
        curr_node = self.sink
        curr_score = self.get_node(len(self.col), len(self.row))

        while curr_node != self.source:
            i = curr_node[0]
            j = curr_node[1]

            u_score = self.get_node(i, j - 1) + self.indel
            l_score = self.get_node(i - 1, j) + self.indel

            if self.col[i - 1] == self.row[j - 1]:
                d_score = self.get_node(i - 1, j - 1) + self.match
            else:
                d_score = self.get_node(i - 1, j - 1) + self.mis

            if d_score == curr_score:
                col_seq.append(self.col[i - 1])
                row_seq.append(self.row[j - 1])
                curr_score = self.get_node(i - 1, j - 1)
                curr_node = (i - 1, j - 1)
                continue

            elif u_score == curr_score:
                col_seq.append('-')
                row_seq.append(self.row[j - 1])
                curr_score = u_score - self.indel
                curr_node = (i, j - 1)
                continue

            elif l_score == curr_score:
                col_seq.append(self.col[i - 1])
                row_seq.append('-')
                curr_score = l_score - self.indel
                curr_node = (i - 1, j)
                continue

        col_seq.reverse()
        row_seq.reverse()

        col_align = ''.join(col_seq)
        row_align = ''.join(row_seq)

        return [self.get_node(len(self.col), len(self.row)),
                len(col_align), col_align, row_align]

    def get_node(self, i, j):
        """
        Retrieve Node at given index
        :param int i: Column index
        :param int j: Row Index
        :return tuple: Node represented as (col, row)
        """
        curr_row = self.matrix[j]
        return curr_row[i]


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


def global_alignment(seq1, seq2, match, mis, indel, to_print):
    """
    Perform Global Alignment and output Score + length
    :param str seq1: DNA Sequence
    :param str seq2: DNA Sequence
    :param int match: Score value for match
    :param int mis: Score penalty for mismatch
    :param int indel: Score penalty for indel
    :param bool to_print: Outputs alignment if true
    """
    scoring_matrix = AlignmentGraph(seq1, seq2, match, mis, indel)
    alignment = scoring_matrix.back_trace()
    print("Score: " + str(alignment[0]))
    print("Length: " + str(alignment[1]))

    if to_print:
        print(alignment[2])
        print(alignment[3])


def main():
    parser = argparse.ArgumentParser(
        description='Calculate optimal global alignment of two seqs')

    parser.add_argument('sequence_file', help='Input sequence file')
    parser.add_argument('-m', '--match', type=int, help='Match value')
    parser.add_argument('-s', '--mismatch', type=int,
                        required=True, help='Mismatch value')
    parser.add_argument('-d', '--indel', type=int,
                        required=True, help='Indel value')
    parser.add_argument('-a', action='store_true', help='Output Alignment')
    args = parser.parse_args()

    match = args.match
    mis = args.mismatch
    indel = args.indel
    to_print = args.a

    seq_list = parse_multiseq(args.sequence_file)
    seq1 = seq_list[0]
    seq2 = seq_list[1]

    global_alignment(seq1, seq2, match, mis, indel, to_print)


if __name__ == '__main__':
    main()
