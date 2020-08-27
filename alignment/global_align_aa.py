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

    Amino Acid scoring through BLOSUM62 substitution
    """
    def __init__(self, col, row, indel):
        """
        Create and score alignment graph
        :param string col: DNA Sequence
        :param string row: DNA Sequence
        :param int indel: Score penalty for indel
        """
        self.matrix = []
        self.col = col
        self.row = row
        self.indel = indel
        self.source = (0, 0)
        self.sink = (len(col), len(row))
        self.aa_index = blosum62()[0]
        self.score_matrix = blosum62()[1]

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

                        scores.append(self.get_node(i - 1, j - 1)
                                      + self.get_score(col[i - 1], row[j - 1]))
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
            d_score = (self.get_node(i - 1, j - 1) +
                       self.get_score(self.col[i - 1], self.row[j - 1]))

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

    def get_score(self, col_aa, row_aa):
        """
        Get score between two amino acids
        :param str col_aa: Amino Acid in column
        :param str row_aa: Amino Acid in row
        :return int: Score
        """
        return self.score_matrix[self.aa_index[col_aa]][self.aa_index[row_aa]]


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


def global_alignment(seq1, seq2, indel, to_print):
    """
    Perform Global Alignment and output Score + length
    :param str seq1: DNA Sequence
    :param str seq2: DNA Sequence
    :param int indel: Score penalty for indel
    :param bool to_print: Outputs alignment if true
    """
    scoring_matrix = AlignmentGraph(seq1, seq2, indel)
    alignment = scoring_matrix.back_trace()
    print("Score: " + str(alignment[0]))
    print("Length: " + str(alignment[1]))

    if to_print:
        print(alignment[2])
        print(alignment[3])


def blosum62():
    """
    Create blosum62 scoring matrix
    :return dict, list: Dict- Index for Amino Acid location in matrix
                        List- List of lists containing difference scores
                              score denoted by list[row_index][col_index]
    """
    bl_str = ("   A  C  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  Y\n"
              "A  4  0 -2 -1 -2  0 -2 -1 -1 -1 -1 -2 -1 -1 -1  1  0  0 -3 -2\n"
              "C  0  9 -3 -4 -2 -3 -3 -1 -3 -1 -1 -3 -3 -3 -3 -1 -1 -1 -2 -2\n"
              "D -2 -3  6  2 -3 -1 -1 -3 -1 -4 -3  1 -1  0 -2  0 -1 -3 -4 -3\n"
              "E -1 -4  2  5 -3 -2  0 -3  1 -3 -2  0 -1  2  0  0 -1 -2 -3 -2\n"
              "F -2 -2 -3 -3  6 -3 -1  0 -3  0  0 -3 -4 -3 -3 -2 -2 -1  1  3\n"
              "G  0 -3 -1 -2 -3  6 -2 -4 -2 -4 -3  0 -2 -2 -2  0 -2 -3 -2 -3\n"
              "H -2 -3 -1  0 -1 -2  8 -3 -1 -3 -2  1 -2  0  0 -1 -2 -3 -2  2\n"
              "I -1 -1 -3 -3  0 -4 -3  4 -3  2  1 -3 -3 -3 -3 -2 -1  3 -3 -1\n"
              "K -1 -3 -1  1 -3 -2 -1 -3  5 -2 -1  0 -1  1  2  0 -1 -2 -3 -2\n"
              "L -1 -1 -4 -3  0 -4 -3  2 -2  4  2 -3 -3 -2 -2 -2 -1  1 -2 -1\n"
              "M -1 -1 -3 -2  0 -3 -2  1 -1  2  5 -2 -2  0 -1 -1 -1  1 -1 -1\n"
              "N -2 -3  1  0 -3  0  1 -3  0 -3 -2  6 -2  0  0  1  0 -3 -4 -2\n"
              "P -1 -3 -1 -1 -4 -2 -2 -3 -1 -3 -2 -2  7 -1 -2 -1 -1 -2 -4 -3\n"
              "Q -1 -3  0  2 -3 -2  0 -3  1 -2  0  0 -1  5  1  0 -1 -2 -2 -1\n"
              "R -1 -3 -2  0 -3 -2  0 -3  2 -2 -1  0 -2  1  5 -1 -1 -3 -3 -2\n"
              "S  1 -1  0  0 -2  0 -1 -2  0 -2 -1  1 -1  0 -1  4  1 -2 -3 -2\n"
              "T  0 -1 -1 -1 -2 -2 -2 -1 -1 -1 -1  0 -1 -1 -1  1  5  0 -2 -2\n"
              "V  0 -1 -3 -2 -1 -3 -3  3 -2  1  1 -3 -2 -2 -3 -2  0  4 -3 -1\n"
              "W -3 -2 -4 -3  1 -2 -2 -3 -3 -2 -1 -4 -4 -2 -3 -3 -2 -3 11  2\n"
              "Y -2 -2 -3 -2  3 -3  2 -1 -2 -1 -1 -2 -3 -1 -2 -2 -2 -1  2  7"
              )

    bl_list = [row.strip() for row in bl_str.split('\n')]
    aa_index = {aa: i for i, aa in enumerate(bl_list[0].split())}

    blosum_matrix = []
    for i in range(1, len(bl_list)):
        curr_row = list(map(int, bl_list[i].split()[1:]))
        blosum_matrix.append(curr_row)

    return aa_index, blosum_matrix


def main():
    parser = argparse.ArgumentParser(
        description='Calculate optimal global alignment of two seqs')

    parser.add_argument('sequence_file', help='Input sequence file')
    parser.add_argument('-d', '--indel', type=int,
                        required=True, help='Indel value')
    parser.add_argument('-a', action='store_true', help='Output Alignment')
    args = parser.parse_args()

    indel = args.indel
    to_print = args.a

    seq_list = parse_multiseq(args.sequence_file)
    seq1 = seq_list[0]
    seq2 = seq_list[1]

    global_alignment(seq1, seq2, indel, to_print)


if __name__ == '__main__':
    main()
