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
