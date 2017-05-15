import random

kmers = []

def allklength(k, alphabet):
    allklength_helper(k, alphabet, "")
    sequence_map = {}
    for i in range(len(kmers)):
        sequence_map[kmers[i]] = i;
    return sequence_map

def allklength_helper(k, alphabet, prefix):
    if k == 0:
        global kmers
        kmers.append(prefix)
        return
    for i in range(len(alphabet)):
        new_prefix = prefix + alphabet[i]
        allklength_helper(k-1, alphabet, new_prefix)

def calculate_distribution(input_sequence, k, alphabet):
    matrix_dim = len(alphabet)^k;
    transition_matrix = [[0*matrix_dim]*matrix_dim]
    sequence_map = allklength(k, alphabet)
    for i in range(len(input_sequence)-k-2):
        current = input_sequence[i:i+k]
        next = input_sequence[i+1:i+k+1]
        transition_matrix[sequence_map[current], sequence_map[next]] += 1
    for i in range(matrix_dim):
        norm_factor = sum(transition_matrix[i, 0:matrix_dim-1])
        for j in range(matrix_dim):
            transition_matrix[i, j] = float[transition_matrix[i, j]]/norm_factor
    return transition_matrix, sequence_map
        
def deletionChannel(input_sequence, deletion_rate):
    for i in range(len(input_sequence)):
        x = random.random()
        if x < deletion_rate:
            input_sequence[i] = 'X'
    noisy = input_sequence.replace('X', "")
    return noisy


def denoiseSequence(input_sequence, k, alphabet, deletion_rate):
    noisy = deletionChannel(input_sequence, deletion_rate)
    pi, seq_map = calculate_distribution(noisy, k, alphabet)
    denoised = noisy[:k]
    for i in range(k-1, len(noisy)-1-k):
        base_prob = 1-deletion_rate
        next_char = noisy[i+1]
        subseq = noisy[i+1-k:i+k+1]
        for j in range(k+1):
            base_prob *= pi[seq_map[subseq[j:j+k]], seq_map[subseq[j+1:j+k+1]]]
        max_prob = base_prob
        for a in alphabet:
            insert_prob = deletion_rate
            new_subseq = subseq[:k+1]+a+subseq[k+1:]
            for m in range(k+2):
                insert_prob *= pi[seq_map[new_subseq[m:m+k]], seq_map[new_subseq[m+1:m+k+1]]]
            if insert_prob > max_prob:
                max_prob = insert_prob
                next_char = a + noisy_char[i+1]
        denoised += next_char
    return input_sequence, noisy, denoised

input = '0101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101'
k = 1
alphabet = ['0', '1']
deletion_rate = 0.2
input, noisy, denoised = denoiseSequence(input, k, alphabet, deletion_rate)
