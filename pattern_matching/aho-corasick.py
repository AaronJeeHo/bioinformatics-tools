"""
Author: Aaron Ho
Python Version: 3.7
"""
import argparse
import time
import csv


class AhoCorasick:
    """
    Class that creates a trie to pattern match
    queries to a database using the aho-corasick algorithm
    """
    def __init__(self):
        """
        Initialize matching trie
        """
        self.nodes = [{'A': None, 'C': None, 'G': None, 'T': None}]
        self.suffix_links = ['']
        self.fail_links = [0]
        self.dict_links = {}

    def insert(self, data):
        """
        Insert database item into trie
        :param str data: sequence from database
        """
        curr_node = 0
        for i in data:

            # add node if doesn't exist
            if self.nodes[curr_node][i] is None:
                self.nodes[curr_node][i] = len(self.nodes)
                self.nodes.append({'A': None, 'C': None, 'G': None, 'T': None})
                self.suffix_links.append(self.suffix_links[curr_node] + i)
                curr_node = len(self.nodes) - 1
            else:
                curr_node = self.nodes[curr_node][i]

        # create dictionary link once finished
        self.dict_links[curr_node] = self.suffix_links[curr_node]

    def set_fail_links(self):
        """
        Sets failure links for trie
        """
        # reset fail links to avoid issues
        self.fail_links = [0]

        for node in range(1, len(self.nodes)):

            suff = self.suffix_links[node][1:]
            while len(self.fail_links) <= node:

                # if longest suffix is empty
                if suff is '':
                    self.fail_links.append(0)
                else:
                    curr_node = 0
                    for i in suff:

                        # check if suffix path exists
                        if self.nodes[curr_node][i] is None:
                            curr_node = None
                            break
                        else:
                            curr_node = self.nodes[curr_node][i]

                    if curr_node is not None:
                        self.fail_links.append(curr_node)
                    else:
                        # shorten suffix by 1
                        suff = suff[1:]


def read_genome(f_name):
    """
    Get genome database from file
    :param str f_name: Database file name
    :return str: Database Genome
    """
    with open(f_name, "r") as file:
        # next(file)
        seq_data = ''.join([line.rstrip() for line in file
                            if not line.startswith('>')])
    return seq_data


def read_dict(f_name):
    """
    Retrieve query sequences
    :param str f_name: Query file name
    :return list: List of seqs to query
    """
    with open(f_name, "r") as file:
        d_lines = [line.rstrip() for line in file]
    return d_lines


def write_m_index(f_name, i_list):
    """
    Create file with query matches and
    their corresponding index
    :param str f_name: Name of output file
    :param list i_list: List containing indeces
    """
    with open(f_name, "w") as file:
        tsv_writer = csv.writer(file, delimiter='\t')
        tsv_writer.writerow(['Matched Strand', 'Start', 'End'])
        for i in i_list:
            tsv_writer.writerow([i[0], str(i[1]), str(i[2])])


def write_estats(f_name, e_matches, num_matches):
    """
    Create file with expected values and
    match stats after pattern matching
    :param str f_name: Name of output file
    :param dict e_matches: Dictionary containing number
                            of expected matches per query
    :param dict num_matches: Dictionary containing actual number
                            of matches per query
    """
    with open(f_name, "w") as file:
        total_expected = sum(e_matches[e] for e in e_matches)
        total_matches = sum(num_matches[m] for m in num_matches)

        tsv_writer = csv.writer(file, delimiter='\t')
        file.write(f'Total Expected Matches: {total_expected}\n')
        file.write(f'Actual Number of Matches: {total_matches}\n\n')

        tsv_writer.writerow(['L-mer', 'Expected Matches', 'Actual Matches'])
        for i in e_matches:
            tsv_writer.writerow([i, str(e_matches[i]), str(num_matches[i])])


def e_values(query_list, genome):
    """
    Create dictionary containing expected
    number of matches per query
    :param list query_list: List of queries
    :param str genome: Genome sequence
    :return dict: Dictionary containing E-values
    """
    expected_matches = {}
    seq_len = len(genome)

    for seq in query_list:
        lmer_len = len(seq)

        if seq not in expected_matches:
            expected_matches[seq] = (seq_len - lmer_len + 1) / (4 ** lmer_len)
        else:
            expected_matches[seq] += (seq_len - lmer_len + 1) / (4 ** lmer_len)

    return expected_matches


def aho_corasick_automaton(query_list):
    """
    Create Aho-Corasick Trie
    :param list query_list: Query to be added into db
    :return AhoCorasick: Aho-Corasick Trie
    """
    aho_trie = AhoCorasick()

    print('Creating Aho-Corasick Automaton')
    for data in query_list:
        aho_trie.insert(data)

    aho_trie.set_fail_links()

    print(f'Created with {len(aho_trie.dict_links)}'
          f' seqs and {len(aho_trie.nodes)} nodes\n')
    return aho_trie


def pattern_matching(aho_trie, genome):
    """
    Match patterns using the Aho-Corasick Trie
    :param AhoCorasick aho_trie: Aho-Corasick Trie
    :param str genome: Genome Sequence
    :return Dict,List : pattern_counts-dict containing number of matches
                        match_index-List containing match locations
    """
    curr_node = start_index = i = 0
    pattern_counts = {}
    match_index = []

    for pattern in aho_trie.dict_links:
        pattern_counts[aho_trie.dict_links[pattern]] = 0

    start = time.perf_counter()
    while start_index < len(genome) and i < len(genome):

        # checking for matching failure
        if aho_trie.nodes[curr_node][genome[i]] is None:
            fail_node = aho_trie.fail_links[curr_node]

            # if failed at root
            if curr_node == 0:
                i += 1
                start_index = i
                curr_node = fail_node
            else:
                start_index = i - len(aho_trie.suffix_links[fail_node])
                curr_node = fail_node
        else:
            curr_node = aho_trie.nodes[curr_node][genome[i]]

            if curr_node in aho_trie.dict_links:
                match_index.append((aho_trie.suffix_links[curr_node], start_index, i))
                pattern_counts[aho_trie.suffix_links[curr_node]] += 1

            i += 1

    end = time.perf_counter()
    print(f'\nMatching completed in {end - start} seconds')
    return pattern_counts, match_index


def main():
    parser = argparse.ArgumentParser(description='Aho-Corasick Matching Algorithm')
    parser.add_argument('-d', '--database', type=str, required=True, help='Database Genome')
    parser.add_argument('-q', '--query', type=str, required=True, help='Query to Match')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output File')
    args = parser.parse_args()

    genome = read_genome(args.database)
    data_list = read_dict(args.query)

    aho_trie = aho_corasick_automaton(data_list)
    num_matches, m_index = pattern_matching(aho_trie, genome)
    expected_matches = e_values(data_list, genome)

    write_m_index(f'{args.output}.tsv', m_index)
    write_estats(f'{args.output}_stats.tsv', expected_matches, num_matches)

    print(f"Matches written to '{args.output}.tsv'")
    print(f"Stats written to '{args.output}_stats.tsv'")


if __name__ == '__main__':
    main()
