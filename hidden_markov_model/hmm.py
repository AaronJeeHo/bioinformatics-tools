"""
Author: Aaron Ho
Python Version: 3.7
"""
import argparse


class HMM:
    """
    Class that creates hidden markov model
    """
    def __init__(self, states, i_prob, t_prob, symbols, e_prob):
        """
        Initialize HMM using parsed and formatted data
        :param list states: List of states
        :param dict i_prob: Initial state probability
                            Key(str) : Initial state
                            Value(int): Initial prob
        :param nested dict t_prob: Transition probability
                                    Key(str): Current state
                                    Value(dict): Dict with transition prob
                                        Key(str): Transition state
                                        Value(int): Transition prob

        :param list symbols: List of possible symbols emitted
        :param nested dict e_prob: Emission probability
                                    Key(str): Current state
                                    Value(dict): Dict with emission prob
                                        Key(str): Symbol
                                        Value(int): Emission prob
        """
        self.states = states
        self.i_prob = i_prob
        self.t_prob = t_prob
        self.symbols = symbols
        self.e_prob = e_prob

    def score_edge(self, k, i, s):
        """
        Returns probability state k transitions to state i
        and emits symbol s
        :param str k: Current state
        :param str i: Next state
        :param str s: Emitted symbol at next state
        :return int: Edge score
        """
        return self.t_prob[k][i] * self.e_prob[i][s]


class HMMAutomaton:
    """
    Class that initializes HMM Automaton for
    dynamic programming based algorithms
    """
    def __init__(self, hmm, emit):
        """
        Initializes a Manhattan style HMM graph
        :param HMM hmm: HMM object
        :param str emit: Emitted sequence
        """
        self.hmm = hmm
        self.graph = {}
        self.source = 1.0
        self.sink = 0.0
        self.emit_list = list(emit)

        for state in hmm.states:
            self.graph[state] = [0.0] * len(emit)
            self.graph[state][0] = (self.source * hmm.i_prob[state]
                                    * hmm.e_prob[state][self.emit_list[0]])

    def viterbi(self):
        """
        Using the viterbi  algorithm,
        determines most probable path for HMM
        to output the given string
        :return str: Most probable path
        """
        for i in range(1, len(self.emit_list)):

            for state in self.hmm.states:
                score_list = []

                for prev_state in self.hmm.states:
                    score_list.append(
                        self.graph[prev_state][i - 1] * self.hmm.score_edge(
                            prev_state, state, self.emit_list[i]))

                self.graph[state][i] = max(score_list)

        max_prob = 0.0
        for f_state in self.hmm.states:
            # curr_prob = self.graph[f_state][len(self.emit_list) - 1]
            if max_prob < self.graph[f_state][len(self.emit_list) - 1]:
                max_prob = self.graph[f_state][len(self.emit_list) - 1]
                self.sink = f_state

        # Backtrack
        path = [self.sink]

        curr_state = self.sink
        curr_score = max_prob
        curr_index = len(self.emit_list) - 1
        while 0 < curr_index:
            for prev_state in self.hmm.states:
                prev_score = self.graph[prev_state][curr_index - 1]
                prev_weight = self.hmm.score_edge(
                    prev_state, curr_state, self.emit_list[curr_index])

                if prev_score * prev_weight == curr_score:
                    curr_state = prev_state
                    curr_score = prev_score

            path.append(curr_state)
            curr_index -= 1

        path.reverse()

        return ''.join(path)


def read_combined(parse_order, hmm_file):
    """
    Parse file containing HMM data
    :param str parse_order: 5 character string denoting parse order
                            'q'- states
                            'i'- initial state probabilities
                            't'- state transition probabilities
                            's'- symbols emitted
                            'e'- emission probabilities

    :param str hmm_file: file containing HMM data separated with '-'
                        FORMAT: whitespace separated values
                            -states: Single line
                                eg. A B
                            -initial prob: Header and data
                                eg. A   B
                                    0.5 0.5
                            -transition prob: Matrix formatted data
                                eg.     A   B
                                    A   0.641   0.359
                                    B   0.729   0.271
                            -symbols: Single line
                            -emission prob: Matrix formatted data

    :return dict: dictionary containing parsed HMM data
                    KEYS: states, i_prob, t_prob, symbols, e_prob
    """
    order = list(parse_order)
    parsed_data = {}
    raw_data = {}
    curr_data = []
    section = []

    with open(hmm_file) as file:
        for line in file:
            if line.startswith('-'):
                section.append(curr_data)
                curr_data = []
            else:
                curr_data.append(line.strip())

        if len(curr_data) < 5:
            section.append(curr_data)

    count = 0
    for i in order:
        raw_data[i] = section[count]
        count += 1

    # process state
    parsed_data['states'] = raw_data['q'][0].split()

    # process i_prob
    i_data = {}
    i_header = raw_data['i'][0].split()
    i_vals = raw_data['i'][1].split()

    i_count = 0
    for i_state in i_header:
        i_data[i_state] = float(i_vals[i_count])
        i_count += 1
    parsed_data['i_prob'] = i_data

    # process t_prob
    t_data = {}
    t_header = raw_data['t'][0].split()

    for t_row in range(1, len(raw_data['t'])):
        curr_row = raw_data['t'][t_row].split()
        curr_state = curr_row[0]
        state_probs = {}

        for s_col in range(len(t_header)):
            state_probs[t_header[s_col]] = float(curr_row[s_col + 1])

        t_data[curr_state] = state_probs
    parsed_data['t_prob'] = t_data

    # process symbols
    parsed_data['symbols'] = raw_data['s'][0].split()

    # process e_prob
    e_data = {}
    e_header = raw_data['e'][0].split()

    for e_row in range(1, len(raw_data['e'])):
        curr_row = raw_data['e'][e_row].split()
        curr_state = curr_row[0]
        emit_probs = {}

        for e_col in range(len(e_header)):
            emit_probs[e_header[e_col]] = float(curr_row[e_col + 1])

        e_data[curr_state] = emit_probs
    parsed_data['e_prob'] = e_data

    return parsed_data


def make_hmm(parsed_data):
    """
    Create HMM after reading and parsing HMM data
    :param dict parsed_data: Parsed hmm data
    :return: HMM object
    """
    return HMM(parsed_data['states'], parsed_data['i_prob'],
               parsed_data['t_prob'], parsed_data['symbols'],
               parsed_data['e_prob'])


def prob_path(hmm, path):
    """
    Probability an hmm follows a given path
    :param HMM hmm: HMM object
    :param str path: Path of state transitions
    :return int: Probability of path
    """
    p_list = list(path)
    prob = 1 * hmm.i_prob[p_list[0]]

    for i in range(1, len(p_list)):
        prob = prob * hmm.t_prob[p_list[i - 1]][p_list[i]]
    return prob


def prob_path_emission(hmm, path, emit):
    """
    Find the probability a given path
    would emit a specific outcome
    :param HMM hmm: HMM object
    :param str path: Path of states
    :param str emit: Specific output
    :return int: Probability of outcome
    """
    p_list = list(path)
    e_list = list(emit)
    prob = 1
    for i in range(len(p_list)):
        prob = prob * hmm.e_prob[p_list[i]][e_list[i]]
    return prob


def main():
    pass


if __name__ == '__main__':
    main()
