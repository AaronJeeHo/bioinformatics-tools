"""
Author: Aaron Ho
Python Version: 3.7
"""
import argparse


class AlignmentGraph:
    """
    Class that creates dynamic programming graph
    and handles alignment scoring through
    Smith-Waterman algorithm.

    Amino Acid scoring through PAM250 substitution
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
        self.max_score = 0
        self.max_node = None
        self.aa_index = pam250()[0]
        self.score_matrix = pam250()[1]

        # fill matrix with empty rows
        for j in range(len(row) + 1):
            self.matrix.append([])

        # start scoring
        for j in range(len(row) + 1):
            curr_row = self.matrix[j]

            if j == 0:
                for i in range(len(col) + 1):
                    curr_row.append(0)
            else:
                for i in range(len(col) + 1):
                    if i == 0:
                        curr_row.append(0)
                    else:
                        scores = []  # [diag , left , up]

                        scores.append(self.get_node(i - 1, j - 1)
                                      + self.get_score(col[i - 1], row[j - 1]))
                        scores.append(self.get_node(i, j - 1) + self.indel)
                        scores.append(self.get_node(i - 1, j) + self.indel)

                        curr_max = max(scores)

                        if curr_max > 0:
                            curr_row.append(curr_max)
                        else:
                            curr_row.append(0)

                        if curr_max > self.max_score:
                            self.max_score = curr_max
                            self.max_node = (i, j)

    def back_trace(self):
        """
        Backtrack and find alignment path
        :return list: List containing [MAX SCORE, ALIGNMENT LEN,
                        COLUMN ALIGNMENT, ROW ALIGNMENT]
        """
        col_seq = []
        row_seq = []
        curr_score = self.max_score
        curr_node = self.max_node

        while curr_score != 0:
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

        return [self.max_score, len(col_align), col_align, row_align]

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


def local_alignment(seq1, seq2, indel, to_print):
    """
    Perform Local Alignment and output Score + length
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


def pam250():
    """
    Create pam250 scoring matrix
    :return dict, list: Dict- Index for Amino Acid location in matrix
                        List- List of lists containing difference scores
                              score denoted by list[row_index][col_index]
    """
    pam_st = ("   A  C  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  Y\n"
              "A  2 -2  0  0 -3  1 -1 -1 -1 -2 -1  0  1  0 -2  1  1  0 -6 -3\n"
              "C -2 12 -5 -5 -4 -3 -3 -2 -5 -6 -5 -4 -3 -5 -4  0 -2 -2 -8  0\n"
              "D  0 -5  4  3 -6  1  1 -2  0 -4 -3  2 -1  2 -1  0  0 -2 -7 -4\n"
              "E  0 -5  3  4 -5  0  1 -2  0 -3 -2  1 -1  2 -1  0  0 -2 -7 -4\n"
              "F -3 -4 -6 -5  9 -5 -2  1 -5  2  0 -3 -5 -5 -4 -3 -3 -1  0  7\n"
              "G  1 -3  1  0 -5  5 -2 -3 -2 -4 -3  0  0 -1 -3  1  0 -1 -7 -5\n"
              "H -1 -3  1  1 -2 -2  6 -2  0 -2 -2  2  0  3  2 -1 -1 -2 -3  0\n"
              "I -1 -2 -2 -2  1 -3 -2  5 -2  2  2 -2 -2 -2 -2 -1  0  4 -5 -1\n"
              "K -1 -5  0  0 -5 -2  0 -2  5 -3  0  1 -1  1  3  0  0 -2 -3 -4\n"
              "L -2 -6 -4 -3  2 -4 -2  2 -3  6  4 -3 -3 -2 -3 -3 -2  2 -2 -1\n"
              "M -1 -5 -3 -2  0 -3 -2  2  0  4  6 -2 -2 -1  0 -2 -1  2 -4 -2\n"
              "N  0 -4  2  1 -3  0  2 -2  1 -3 -2  2  0  1  0  1  0 -2 -4 -2\n"
              "P  1 -3 -1 -1 -5  0  0 -2 -1 -3 -2  0  6  0  0  1  0 -1 -6 -5\n"
              "Q  0 -5  2  2 -5 -1  3 -2  1 -2 -1  1  0  4  1 -1 -1 -2 -5 -4\n"
              "R -2 -4 -1 -1 -4 -3  2 -2  3 -3  0  0  0  1  6  0 -1 -2  2 -4\n"
              "S  1  0  0  0 -3  1 -1 -1  0 -3 -2  1  1 -1  0  2  1 -1 -2 -3\n"
              "T  1 -2  0  0 -3  0 -1  0  0 -2 -1  0  0 -1 -1  1  3  0 -5 -3\n"
              "V  0 -2 -2 -2 -1 -1 -2  4 -2  2  2 -2 -1 -2 -2 -1  0  4 -6 -2\n"
              "W -6 -8 -7 -7  0 -7 -3 -5 -3 -2 -4 -4 -6 -5  2 -2 -5 -6 17  0\n"
              "Y -3  0 -4 -4  7 -5  0 -1 -4 -1 -2 -2 -5 -4 -4 -3 -3 -2  0 10"
              )

    pam_list = [row.strip() for row in pam_st.split('\n')]
    aa_index = {aa: i for i, aa in enumerate(pam_list[0].split())}

    pam_matrix = []
    for i in range(1, len(pam_list)):
        curr_row = list(map(int, pam_list[i].split()[1:]))
        pam_matrix.append(curr_row)

    return aa_index, pam_matrix


def main():
    parser = argparse.ArgumentParser(
        description='Calculate optimal local alignment of two seqs')
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

    local_alignment(seq1, seq2, indel, to_print)


if __name__ == '__main__':
    main()
