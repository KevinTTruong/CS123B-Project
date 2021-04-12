gap_char = '-'
filename = ""
instructions = "Instructions:\n" \
                "* = mutation\n- = gap\n" \
                "-The top line indicates gaps/mutation for the OVERALL alignment, while the gaps/mutation\n" \
                "between sequences for each PAIR.\n" \
                "-For the \"Potential mutated\" section, the \"# seq\" indicates how many sequences contain\n" \
                "each INDIVIDUAL function (listed just below) and the #s following indicate possible mutation\n" \
                "locations (the numbers header should help for reference).\n" \
                "-For example, \"[func] was found in [X] sequences and could either be the [nucleotide] at\n" \
                "[location] or [another nucleotide] at [another location]\".\n" \

# List of sequences
sequences = []

# key = functionality; value = # of seqs containing func (frequency)
mut = {}

# key = # of seqs containing func (frequency); value = func (inverse of mut)
mut_by_freq = {}

# key = frequency; value = possible mutation locations for frequency
loc_by_freq = {}


# Needleman-Wunsch algorithm
# Match: +1
# Mismatch: -2
# Gap: -2
def align(seq, index1, index2, alt_priority = False):
    matrix = [[0 for x in range(len(seq[index2])+1)] for y in range(len(seq[index1])+1)]
    # Initialize first row/column
    matrix[0][0] = 0
    for i in range(1, len(seq[index1])+1):
        matrix[i][0] = -2 * i
    for j in range(1, len(seq[index2])+1):
        matrix[0][j] = -2 * j
    # Fill matrix
    for x in range(1, max(len(seq[index1])+1, len(seq[index2])+1)):
        if len(seq[index1]) > x-1 and len(seq[index2]) > x-1:
            for i in range(x, len(seq[index1])+1):
                calc = -2
                if seq[index1][i-1] == seq[index2][x-1]:
                    calc = 1
                matrix[i][x] = max(matrix[i-1][x]-2, matrix[i][x-1]-2, matrix[i-1][x-1]+calc)
            for j in range(x, len(seq[index2])+1):
                calc = -2
                if seq[index1][x-1] == seq[index2][j-1]:
                    calc = 1
                matrix[x][j] = max(matrix[x-1][j]-2, matrix[x][j-1]-2, matrix[x-1][j-1]+calc)

    # Traverse and locate gaps
    i = 0
    j = 0
    i_gaps = []
    j_gaps = []
    for x in range(max(len(seq[index1]), len(seq[index2]))):
        if i+1 >= len(seq[index1]) and j+1 < len(seq[index2]):
            i_gaps.append(x+1)  # Recall: gap in matrix for 0
            j += 1
        elif j+1 >= len(seq[index2]) and i+1 < len(seq[index1]):
            j_gaps.append(x+1)
            i += 1
        else:
            k = max(matrix[i+1][j+1], matrix[i][j+1], matrix[i+1][j])
            if not alt_priority:
                if k == matrix[i+1][j+1]:
                    i += 1
                    j += 1
                # Priority route 1
                elif k == matrix[i+1][j]:
                    j_gaps.append(x+1)
                    i += 1
                elif k == matrix[i][j+1]:
                    i_gaps.append(x+1)   # Recall: gap in matrix for 0
                    j += 1
            # Priority route 2
            else:
                if k == matrix[i+1][j+1]:
                    i += 1
                    j += 1
                elif k == matrix[i][j+1]:
                    i_gaps.append(x+1)   # Recall: gap in matrix for 0
                    j += 1
                elif k == matrix[i+1][j]:
                    j_gaps.append(x+1)
                    i += 1

    # Insert gaps into sequences and updates previous ones if needed
    length = range(max(len(seq[index1]), len(seq[index2])))
    if len(i_gaps) > 0:
        i_gap_ptr = len(i_gaps) - 1
        ptr = index1
        original_len = len(seq[ptr])
        for i in reversed(length):
            if i_gaps[i_gap_ptr] > original_len:
                seq[ptr] += gap_char
                i_gap_ptr -= 1
                if i_gap_ptr == -1:
                    break
            elif i == i_gaps[i_gap_ptr]:
                seq[ptr] = seq[ptr][:i] + gap_char + seq[ptr][i:]
                # Exclusive to i: insert gaps to previous sequences
                if ptr > 0:
                    for p in range(ptr):
                        seq[p] = seq[p][:i] + gap_char + seq[p][i:]
                # -----
                i_gap_ptr -= 1
                if i_gap_ptr == -1:
                    break
    # Same algorithm as above, but for the 2nd/latest sequence
    if len(j_gaps) > 0:
        j_gap_ptr = len(j_gaps) - 1
        ptr = index2
        original_len = len(seq[ptr])
        for j in reversed(length):
            if j_gaps[j_gap_ptr] > original_len:
                seq[ptr] += gap_char
                j_gap_ptr -= 1
                if j_gap_ptr == -1:
                    break
            elif j == j_gaps[j_gap_ptr]:
                seq[ptr] = seq[ptr][:j] + gap_char + seq[ptr][j:]
                j_gap_ptr -= 1
                if j_gap_ptr == -1:
                    break
    # Fills any leading gaps to match sequence lengths
    return fill_remaining_gaps(seq, index1, index2)


# Print out mutation locations heading for the overall alignment (the ---***-*-- thing)
def mutation_pointer():
    output = ""
    for x in range(len(seq[0])):
        clean = True
        for y in range(1, len(seq)):
            if seq[0][x] != seq[y][x]:
                output += "*"
                clean = False
                break
        if clean:
            output += '-'
    return output


# Print mutation locations heading for a selected pair of sequences (the ---***-*-- thing)
def pairwise_mutation_pointer(index1, index2):
    output = ""
    for x in range(len(seq[index1])):
        if seq[index1][x] != seq[index2][x]:
            output += "*"
        else:
            output += '-'
    return output


# Calls align() for all sequences in 'seq'
def multi_align(seq, alt_priority = False):
    if len(seq) == 1:
        print("Only one sequence available")
    elif len(seq) == 0:
        print("No sequences available")
    else:
        for x in range(1, len(seq)):
            seq = align(seq, x-1, x, alt_priority)
    return seq


# Fills any leading gaps to match sequence lengths
def fill_remaining_gaps(seq, index1, index2):
    while len(seq[index1]) > len(seq[index2]):
        seq[index2] += gap_char
    while len(seq[index1]) < len(seq[index2]):
        seq[index1] += gap_char
    return seq


# "Main" starts here
#sequences.append("GAGCAGCTGAACAAGCTGATGACCACCCTCCATAGCACCGCACCCCATTTTGTCCGCTGTATTATCCCCAATGAGTTTAAGCAATCGG")
#sequences.append("GAGCAGCTGAACAAGCTGATGACCACCCTCCATAGCACCGCACCCCATTTTGTCCGCTGTATTATCCCCAATGAGTTTAAGCAATCGG")
#sequences.append("GAGCAGCTGAACAAGCTGATGACCACCCTCCACAGCACTGCACCCCATTTTGTCCGCTGTATTGTGCCCAATGAGTTTAAGCAGTCAG")
#sequences.append("GAGCAGCTGAACAAGCTGATGACCACCCTCCATAGCACCGCACCCCATTTTGTCCGCTGTATTATCCCCAATGAGTTTAAGCAATCGG")
#sequences.append("GAGCAGCTGAACAAGCTGATGACCACCCTCCACAGCACTGCACCCCATTTTGTCCGCTGTATTGTGCCCAATGAGTTTAAGCAGTCAG")
#sequences.append("GAGCAGCTGAACAAGCTGATGACCACCCTCCATAGCCGCACCCCATTTTGTCCGCTGTATTATCCCCAATGAGTTTAAGCAATCGG")
#sequences.append("CGCAC")
#sequences.append("GGCTTC")

#sequences = multi_align(sequences)
#for s in sequences:
#    print(s)

# Test cases
# print(mutation_pointer())
# print(print_results())
