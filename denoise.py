import random
import numpy as np

kmers = []

def allklength(k, alphabet):
    allklengthHelper(k, alphabet, "")

def allklengthHelper(k, alphabet, prefix):
    if k == 0:
        global kmers
        kmers.append(prefix)
        return
    for i in range(len(alphabet)):
        new_prefix = prefix + alphabet[i]
        allklengthHelper(k-1, alphabet, new_prefix)

def calculateMatrix(input_sequence, k, alphabet):
    matrix_dim = len(alphabet)**k;
    transition_matrix = np.zeros((matrix_dim, matrix_dim))
    allklength(k, alphabet)
    sequence_map = {}
    for i in range(len(kmers)):
        sequence_map[kmers[i]] = i
    print(sequence_map)
    for i in range(len(input_sequence)-k-2):
        current = input_sequence[i:i+k]
        next = input_sequence[i+1:i+k+1]
        transition_matrix[sequence_map[current], sequence_map[next]] += 1
    for i in range(matrix_dim):
        norm_factor = sum(transition_matrix[i, 0:matrix_dim-1])
        for j in range(matrix_dim):
            transition_matrix[i, j] = transition_matrix[i, j]/(norm_factor+0.0)
    return transition_matrix, sequence_map

def calculateDistribution(input_sequence, k, alphabet):
    context_hist = {}
    for i in range(k, len(input_sequence)-k):
        context = input_sequence[i-k:i]+input_sequence[i+1:i+k+1]
        if context in context_hist:
            context_hist[context] += input_sequence[i]
        else:
            context_hist[context] = [input_sequence[i]]
    return context_hist
        
def erasureChannel(input_sequence, deletion_rate):
    noisy = input_sequence
    for i in range(len(input_sequence)):
        x = random.random()
        if x < deletion_rate:
            noisy = noisy[:i]+'E'+noisy[i+1:]
    return noisy

def deletionChannel(input_sequence, deletion_rate):
    indices = []
    for i in range(len(input_sequence)):
        x = random.random()
        if x < deletion_rate:
            indices += [i]
    noisy = "".join([char for index, char in enumerate(input_sequence) if index not in indices])
    return noisy

def erasureDenoise(input_sequence, k, alphabet, deletion_rate):
    noisy = erasureChannel(input_sequence, deletion_rate)
    contexts = calculateDistribution(noisy, k, alphabet)
    ml = alphabet[0]
    pct = 0
    for l in alphabet:
        if sum([int(x == l) for x in noisy])/len(noisy) > pct:
            pct = sum(noisy == l)/len(noisy)
            ml = l
    most_common = ml
    for j in range(k):
        if noisy[j] == 'E':
            noisy = noisy[:j]+ml+noisy[j+1:]
    for m in range(len(noisy)-k, len(noisy)):
        if noisy[m] == 'E':
            noisy = noisy[:m]+ml+noisy[m+1:]
    erasure_corrections = []
    for i in range(k, len(noisy)-k):
        if noisy[i] == 'E':
            p = 0
            ml = most_common
            context = noisy[i-k:i]+noisy[i+1:i+k+1]
            context_hist = contexts[context]
            for a in alphabet:
                new_p = sum([int(x == a) for x in context_hist])/len(context_hist)
                if new_p > p:
                    p = new_p
                    ml = a
            erasure_corrections += ml
    j = 0
    for i in range(k, len(noisy)-k):
        if noisy[i] == 'E':
            noisy = noisy[:i]+erasure_corrections[j]+noisy[i+1:]
            j += 1
    return noisy

def denoiseSequence(input_sequence, k, alphabet, deletion_rate):
    noisy = deletionChannel(input_sequence, deletion_rate)
    pi, seq_map = calculateMatrix(noisy, k, alphabet)
    denoised = noisy[:k]
    for i in range(k-1, len(noisy)-1-k):
        base_prob = 1-deletion_rate
        next_char = noisy[i+1]
        subseq = noisy[i+1-k:i+k+1]
        for j in range(k):
            base_prob *= pi[seq_map[subseq[j:j+k]], seq_map[subseq[j+1:j+k+1]]]
        max_prob = base_prob
        for a in alphabet:
            insert_prob = deletion_rate
            new_subseq = subseq[:k+1]+a+subseq[k+1:]
            for m in range(k+1):
                insert_prob *= pi[seq_map[new_subseq[m:m+k]], seq_map[new_subseq[m+1:m+k+1]]]
            if insert_prob > max_prob:
                max_prob = insert_prob
                next_char = a + noisy[i+1]
        denoised += next_char
    for i in range(len(noisy)-1-k, len(noisy)):
        denoised += noisy[i]
    return input_sequence, noisy, denoised


n = 100000
rho = 0.01
k = 2
alphabet = ['a', 'c', 't', 'g']
input = open('humangenome.fasta').read()[:n]
erasure_denoised = erasureDenoise(input, k, alphabet, rho) 
print(1-sum([int(input[i] == erasure_denoised[i]) for i in range(len(input))])/(len(input)+0.0))
#input, noisy, denoised = denoiseSequence(input, k, alphabet, rho)

