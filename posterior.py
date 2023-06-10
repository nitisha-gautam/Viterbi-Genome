import sys
from math import log, exp
from compsci260lib import get_fasta_dict


def run_posterior():
    input_file = "bacterial.genome.fasta"
    hmm_file = "HMM.methanococcus.txt"

    posterior = posterior_decoding(input_file, hmm_file)

    # Report the first and last ten segments in your decoding
    #
    print("The last 10 segments are: ")
    for x in posterior[:10]:
        print(x)
    print("")
    print("The first 10 segments are: ")
    for x in posterior[-10:]:
        print(x)
    print("")

    countState = {}

    for x in posterior:
        if x['state'] in countState.keys():
            countState[x['state']] += 1
        else:
            countState[x['state']] = 1

    print(countState)
    #


def posterior_decoding(input_file, hmm_file):
    """
    Calculate the posterior decoding and return the decoded segments.

    input_file (str): path to input fasta file
    hmm_file (str): path to HMM file

    Returns:
        A list of dictionaries of segments in each state (1-indexed).
        An example output may look like:

        [
            {â€˜startâ€™: 1, â€˜endâ€™: 12, â€˜stateâ€™: â€˜state2â€™},
            {â€˜startâ€™: 13, â€˜endâ€™: 20, â€˜stateâ€™: â€˜state1â€™},
            ...
        ]

    """

    # Open HMM description file
    try:
        f_hmm_file = open(hmm_file, 'r')
    except IOError:
        print("IOError: Unable to open HMM file: %s." % hmm_file)
        print("Exiting.")
        sys.exit(1)

    # Read the state names
    states = f_hmm_file.readline().split()
    K = len(states)

    # Read the initial probability vector (and log transform entries)
    probs = f_hmm_file.readline().split()
    initial_probs = [log(float(prob)) for prob in probs]

    # Read the transition matrix (and log transform entries)
    transitions = [None for _ in range(K)]
    for i in range(K):
        matrix_row_arry = f_hmm_file.readline().split()
        matrix_row = [log(float(trans_prob)) for trans_prob in matrix_row_arry]
        transitions[i] = matrix_row

    # Read the emitted symbols
    emitted_symbols_list = f_hmm_file.readline().split()
    emitted_symbols = {symbol: index for index, symbol in enumerate(emitted_symbols_list)}

    # Read the emission probability matrix (and log transform entries)
    emit_probs = [None for _ in range(K)]
    for i in range(K):
        matrix_row_arry = f_hmm_file.readline().split()
        matrix_row = [log(float(emit_prob)) for emit_prob in matrix_row_arry]
        emit_probs[i] = matrix_row

    f_hmm_file.close()

    seq_dict = get_fasta_dict(input_file)
    emit_str = list(seq_dict.values())[0]  # there's only 1

    print("Done reading sequence of length ", len(emit_str))

    # Run the forward algorithm
    forward = run_forward(states, initial_probs, transitions, emitted_symbols,
                          emit_probs, emit_str)

    # Run the backward algorithm
    backward = run_backward(states, initial_probs, transitions,
                            emitted_symbols, emit_probs, emit_str)

    # Calculate the posterior probabilities
    # Initializing the posterior 2D matrices
    posterior = [[float(0) for _ in range(K)] for _ in range(len(emit_str))]
    for i in range(len(emit_str)):
        # Did not normalize the probabilities (i.e., did not divide by P(X)),
        # because we will only use these probabilities to compare
        # posterior[i][0] versus posterior[i][1]
        for k in range(K):
            posterior[i][k] = forward[i][k] + backward[i][k]

    # Create the list of decoded segments to return
    #
    str = posterior[len(emit_str) - 1]
    start = max(str)

    startPosition = str.index(start)
    end = len(emit_str)
    i = len(emit_str) - 1

    # trackback through the matrix based on the position from the initial highest probaiblity to get the states that
    # have the highest probaiblity of occuring through the length of the sequence
    allPos = []
    while i >= 0:
        start = max(posterior[i - 1])
        nextVal = posterior[i - 1].index(start)
        if nextVal == startPosition:
            startPosition = nextVal
            i -= 1
        else:
            allPos.append({'start': i + 1, 'end': end, 'state': states[startPosition]})
            startPosition = nextVal
            end = i
            i -= 1
        if i == 0:
            allPos.append({'start': i + 1, 'end': end, 'state': states[startPosition]})

    if allPos[len(allPos) - 1]['start'] == allPos[len(allPos) - 2]['start']:
        allPos.pop(len(allPos) - 1)

    return allPos


# this function is used to add logs based on the mathematical formula
def addingLogs(values):
    sum = values[0]
    for x in range(1, len(values)):
        sum = sum + log(1 + exp(values[x] - sum))
    return sum


def run_forward(states, initial_probs, transitions, emitted_symbols,
                emit_probs, emit_str):
    """Calculates the forward (log) probability matrix.

    Arguments:
        states (list of str): list of states as strings
        initial_probs (list of float): list of log(initial probabilities) for each
            state
        transitions (list of list of float): matrix of log(transition probabilities)
        emitted_symbols (dict {str: int}): dictionary mapping emitted symbols to their index
        emit_probs (list of list of float): matrix of log(emission probabilities)
            for each state and emission symbol
        emit_str (str):

    Returns:
        (list of list of floats): matrix of forward (log) probabilities
    """

    K = len(states)
    N = len(emit_str)

    forward = [[float(0) for _ in range(K)] for _ in range(N)]

    # Initialize
    emit_index = emitted_symbols[emit_str[0].upper()]
    for k in range(K):
        forward[0][k] = initial_probs[k] + emit_probs[k][emit_index]

    # Iterate
    previous_stat_probabilities = {}
    for i in range(1, N):
        emit_index = emitted_symbols[emit_str[i].upper()]

        for j in range(K):
            # Compute the forward probabilities for the states
            #
            emission_score = emit_probs[j][emit_index]

            # go through all the possible state values and add the probabilities to a list
            for values in range(K):
                transition = transitions[values][j]  # get transition from previous state to current one
                forward_score = forward[i - 1][values]  # calculate previous forward score
                current_probability = forward_score + transition  # get total probability
                previous_stat_probabilities[values] = current_probability  # add to list

            # find the sum of the probabilities of all states and add using helper method
            sum_prob = addingLogs(list(previous_stat_probabilities.values()))

            # append the sum of the probability sum with the emission score to the matrix
            forward[i][j] = max((emission_score + sum_prob), float('-inf'))
            previous_stat_probabilities.clear()

        #
    return forward


def run_backward(states, initial_probs, transitions, emitted_symbols,
                 emit_probs, emit_str):
    """Calculates the backward (log) probability matrix.

        Arguments:
            states (list of str): list of states as strings
            initial_probs (list of float): list of log(initial probabilities) for
                each state
            transitions (list of list of float): matrix of log(transition
                probabilities)
            emitted_symbols (dict {str: int}): dictionary mapping emitted symbols to their index
            emit_probs (list of list of float): matrix of log(emission
                probabilities) for each state and emission symbol
            emit_str (str):

        Returns:
            (list of list of floats): matrix of backward (log) probabilities
    """

    K = len(states)
    N = len(emit_str)

    backward = [[float(0) for _ in range(K)] for _ in range(N)]

    # Initialize
    for k in range(K):
        backward[N - 1][k] = log(1)  # which is zero, but just to be explicit...

    # Iterate
    previous_stat_probabilities = {}
    for i in range(N - 2, -1, -1):
        emit_index = emitted_symbols[emit_str[i + 1].upper()]

        for j in range(K):
            # Compute the backward probabilities for the states

            # go through each probability for each state
            for values in range(K):
                transition = transitions[j][values]  # find transition probability from last state to current one
                emission_score = emit_probs[values][emit_index]  # get the emission score based on current state
                backward_score = backward[i + 1][values]  # get previous backward score
                current_probability = backward_score + transition + emission_score  # get total probability

                previous_stat_probabilities[values] = current_probability  # append probability of this state to list

            # sum the probabilities of all the possible states
            sum_prob = addingLogs(list(previous_stat_probabilities.values()))

            # add this summation to the matrix
            backward[i][j] = sum_prob
            previous_stat_probabilities.clear()

    return backward


if __name__ == '__main__':
    """Call run_posterior(), do not modify"""
    run_posterior()