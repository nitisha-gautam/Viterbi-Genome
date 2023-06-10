import sys
from math import log
from compsci260lib import get_fasta_dict


def run_viterbi():
    hmm_file = 'HMM.methanococcus.txt'
    input_file = 'bacterial.genome.fasta'
    print("Decoding the Viterbi path for %s using %s" % (input_file, hmm_file))

    vit_ret = viterbi_decoding(input_file, hmm_file)

    # Collect the segment counts for each state
    counts = count_segments(vit_ret)

    # Report the first 10 and last 10 segments of your decoding
    #
    print("The last 10 segments are: ")
    for x in vit_ret[:10]:
        print(x)
    print("")
    print("The first 10 segments are: ")
    for x in vit_ret[-10:]:
        print(x)
    print("")
    #

    # Then report the number of segments that exist in each state.
    #

    # find the total count of segments in each of the states, the total state values in each state over the course of
    # the string, and the length of each of the segments under each state
    totalCounts = {}
    averageLength = {}
    for x in counts:
        print("The " + str(x) + " has " + str(counts[x]) + " number of segments")

    for x in vit_ret:
        if x['state'] in totalCounts:
            totalCounts[x['state']] += ((x['end'] - x['start']) + 1)
            averageLength[x['state']].append((x['end'] - x['start']) + 1)
        else:
            totalCounts[x['state']] = ((x['end'] - x['start']) + 1)
            averageLength[x['state']] = [((x['end'] - x['start']) + 1)]

    # print statistics about the states
    print("Number of nucleotides in each state: ")
    print(totalCounts)

    print("")
    seq_dict = get_fasta_dict(input_file)
    emit_str = list(seq_dict.values())[0]  # there's only 1
    print("Average time spent in each state: ")
    for x in totalCounts:
        print(x + ": " + str(totalCounts[x]/len(emit_str)))

    print("Average segment length for each state: ")
    for x in averageLength:
        print(x + ": " + str(sum(averageLength[x])/len((averageLength[x]))))
    #


def viterbi_decoding(input_file, hmm_file):
    """Calculate the Viterbi decoding of an input sequence

    Arguments:
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

    print("Read a sequence of length", len(emit_str))

    # Create Viterbi table and traceback table
    viterbi = [[0 for _ in range(len(emit_str))] for _ in range(K)]
    pointers = [[0 for _ in range(len(emit_str))] for _ in range(K)]

    # Initialize the first column of the matrix
    for i in range(K):
        in_index = emitted_symbols[emit_str[0].upper()]
        viterbi[i][0] = emit_probs[i][in_index] + initial_probs[i]

    # Build the matrix column by column
    previous_stat_probabilities = {}
    for j in range(1, len(emit_str)):
        in_index = emitted_symbols[emit_str[j].upper()]

        for i in range(K):
            # Compute the entries viterbi[i][j] and pointers[i][j]
            # Tip: Use float('-inf') for the value of negative infinity
            #
            emission_score = emit_probs[i][in_index]

            # go through all the states
            for values in range(K):
                transition = transitions[values][i] # get the transition value from the last state to the current state
                # (values)
                viterbi_score = viterbi[values][j - 1] # find the previous viterbi score
                current_probability = viterbi_score + transition # calculate the overall probability by adding the two
                # scores together

                # add the probability to the list of all possible probabilities for the states
                previous_stat_probabilities[values] = current_probability

            # find the max probability out of the possible probabilities for the states
            max_prob = max(previous_stat_probabilities.values())

            # get the index, or the state, for the max value
            for keys in previous_stat_probabilities.keys():
                if previous_stat_probabilities[keys] == max_prob:
                    indval = keys

            # update the matrix
            pointers[i][j] = indval
            viterbi[i][j] = emission_score + max_prob
            previous_stat_probabilities.clear()
            # 
            pass  # placeholder line
    # Traceback, stored as a list of segments in each state (represented using dictionaries)
    lastValues = [row[len(emit_str) - 1] for row in viterbi]

    # get the start value, or the highest probability from the last column of the matrix
    start = max(lastValues)
    startPosition = lastValues.index(start)
    end = len(emit_str)
    i = len(emit_str) - 1

    # trackback through the matrix based on the position from the initial highest probaiblity to get the states that
    # have the highest probaiblity of occuring through the length of the sequence
    allPos = []
    while i >= 0:
        next = pointers[startPosition][i]
        if next == startPosition:
            startPosition = next
            i -= 1
        else:
            allPos.append({'start': i+1, 'end': end, 'state': states[startPosition]})
            startPosition = pointers[startPosition][i]
            end = i
            i -= 1
        if i == 0:
            allPos.append({'start': i + 1, 'end': end, 'state': states[startPosition]})

    return allPos


def count_segments(vit_ret):
    """Calculate the number of segments appearing in each state of
    the viterbi path

    Arguments:
        vit_ret (list of dicts): dictionary of segments in each state.
            see: return value of viterbi_decoding

    Returns:
        a dictionary mapping states to number of segments of that state. 
        e.g. {'state1': 10, 'state2': 9}
    """

    #

    # count the number of segments found for each state by adding to a dictionary of state names
    segmentCount = {}
    for dicts in vit_ret:
        if dicts['state'] in segmentCount.keys():
            segmentCount[dicts['state']] += 1
        else:
            segmentCount[dicts['state']] = 1

    #
    return segmentCount


if __name__ == '__main__':
    """Call run_viterbi(), do not modify"""
    run_viterbi()