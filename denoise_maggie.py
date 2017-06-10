import random
import numpy as np
import math

global kmers
kmers = []
read_length = 100
num_reads = 1000000

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

def takeReads(genome):
    concatenated = ''
    for i in range(num_reads):
        index = int(random.random()*(len(genome)-read_length))
        concatenated += genome[index:index+read_length]
    return concatenated

def kEntropy(input, k, alphabet):
    return 0

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
    erasure_corrections = []
    for j in range(k):
        if noisy[j] == 'E':
            erasure_corrections += most_common
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
    for m in range(len(noisy)-k, len(noisy)):
        if noisy[m] == 'E':
            erasure_corrections += most_common
    j = 0
    for i in range(len(noisy)):
        if noisy[i] == 'E':
            noisy = noisy[:i]+erasure_corrections[j]+noisy[i+1:]
            j += 1
    return noisy

def levenshtein(source, target):
    if len(source) < len(target):
        return levenshtein(target, source)
    if len(target) == 0:
        return len(source)
    source = np.array(tuple(source))
    target = np.array(tuple(target))
    previous_row = np.arange(target.size + 1)
    for s in source:
        current_row = previous_row + 1
        current_row[1:] = np.minimum(
                current_row[1:],
                np.add(previous_row[:-1], target != s))
        current_row[1:] = np.minimum(
                current_row[1:],
                current_row[0:-1] + 1)
        previous_row = current_row
    return previous_row[-1]

def error(source, target):
    if len(source) < len(target):
        return error(target, source)
    if len(target) == 0:
        return len(source)
    errors = 0
    j = 0
    for i in range(len(source)):
        if j >= len(target):
            errors += 1
        elif source[i] != target[j]:
            errors += 1
        else:
            j += 1
    return errors

def denoiseSequence1(input_sequence, k, alphabet, deletion_rate):
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

def denoiseSequence2(noisy, k, alphabet, deletion_rate):
    #noisy = deletionChannel(input_sequence, deletion_rate)
    adjust = 1
    max_del = 1
    context_hists = []
    for j in range(max_del):
        context_del_hist = {}
        for i in range(k+j, len(noisy)-k):
            context_del = noisy[i-k-j:i-j]+noisy[i:i+k]
            deleted = noisy[i-j:i]
            if context_del in context_del_hist:
                context_del_hist[context_del].append(deleted)
            else:
                context_del_hist[context_del] = [deleted]
        context_hists.append(context_del_hist)
    for j in range(1, max_del):
        allklength(j, alphabet)
        for i in range(k, len(noisy)-k):
            context = noisy[i-k:i+k]
            context_del_hist = context_hists[j]
            context_hist = context_hists[0]
            if context in context_hist and context in context_del_hist:
                if adjust*(deletion_rate**j)*len(context_hist[context]) >= (1.0/adjust)*(1-deletion_rate)*len(context_del_hist[context]):
                    ml = kmers[0]
                    p = 0
                    for a in kmers:
                        new_p = sum([int(x == a) for x in context_hist[context]])/(len(context_hist[context])+0.0)
                        if new_p > p:
                            p = new_p
                            ml = a
                    noisy = noisy[:i]+ml+noisy[i:]
        global kmers
        kmers = []
    return noisy
    
def optimalDenoise(noisy, k, alphabet, rho, alpha):
    for i in range(len(noisy)-1):
        if noisy[i] == noisy[i+1] and rho*alpha**2/(1-alpha) > 1:
            if noisy[i] == alphabet[0]:
                noisy = noisy[:i+1]+alphabet[1]+noisy[i+1:]
            else:
                noisy = noisy[:i+1]+alphabet[0]+noisy[i+1:]
    return noisy
        

def denoiseSequence3(noisy, k, alphabet, rho, l=-1):
    #noisy = deletionChannel(input, rho)
    adjust = 1
    if l == -1:
        l = k
    context_left = {} 
    context_right = {}
    for i in range(k, len(noisy)-k+1):
        if i < len(noisy)-k:
            lcontext = noisy[i-k:i]
            rcontext = noisy[i+1:i+k+1]
            if lcontext in context_left:
                context_left[lcontext] += noisy[i]
            else:
                context_left[lcontext] = [noisy[i]]
            if rcontext in context_right:
                context_right[rcontext] += noisy[i]
            else:
                context_right[rcontext] = [noisy[i]]
    for i in range(k, len(noisy)-k):
        context1 = noisy[i-k+1:i+k+1]
        for a in alphabet:
            context2 = noisy[i-k:i+1]+a+noisy[i+1:i+k+1]
            if a+rcontext[:-1] in context_right and lcontext[1:]+a in context_left:
                #rightp1 = sum([int(x == noisy[i]) for x in context_right[rcontext]])/(len(context_right[rcontext])+0.0)
                #rightp2 = sum([int(x == noisy[i]) for x in context_right[a+rcontext[:-1]]])/(len(context_right[a+rcontext[:-1]])+0.0)
                #leftp1 = sum([int(x == noisy[i+1]) for x in context_left[lcontext]])/(len(context_left[lcontext])+0.0)
                #leftp2 = sum([int(x == noisy[i+1]) for x in context_left[lcontext[1:]+a]])/(len(context_left[lcontext])+0.0)
                if adjust*rho*rightp2 >= (1-rho)*(1.0/adjust)*rightp1 and adjust*rho*leftp2 >= (1-rho)*(1.0/adjust)*leftp1:
                    noisy = noisy[:i+1]+a+noisy[i+1:]
    return noisy

def textDenoise(filename):
    n = 10000
    ks = open(filename, 'r').read()[:n]
    k = int(0.5*math.log(n, 3))
    eps = 0.1
    alphabet = list('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ ')
    noisy = deletionChannel(ks, eps)
    est1 = denoiseSequence2(noisy, k, alphabet, eps)
    f = open(filename+'_denoised_1', 'w')
    f.write(est1)
    f.close()
    f = open(filename+'_noisy', 'w')
    f.write(noisy)
    f.close()

def markovSourceDenoise(a, eps):
    n = 100000
    display = 50
    alphabet = ['+', '-']
    p = random.random()
    k = int(0.5*math.log(n, 3))
    x = ''
    for i in range(n):
        if i == 0:
            if p < 0.5:
                x += '+'
            else:
                x += '-'
        else:
            if p < a:
                if x[i-1] == '+':
                    x += '-'
                else:
                    x += '+'
            else:
                x += x[i-1]
            p = random.random()
    noisy = deletionChannel(x, eps)
    #est1 = denoiseSequence2(noisy, k, alphabet, eps)
    #est2 = denoiseSequence2(noisy, k, alphabet, eps)
    est = optimalDenoise(noisy, k, alphabet, eps, a)
    print 'Setting: alpha = ', a, ', epsilon = ', eps 
    print 'Original: ', x[:display], '(length ', len(x), ' error ', error(x, x)/(n+0.0), ')'
    print 'Noisy: ', noisy[:display], '(length ', len(noisy), ' error ', error(noisy, x)/(n+0.0), ')'
    print 'Denoised: ', est[:display], '(length ', len(est), ' error ', error(est, x)/(n+0.0), ' )'
    #print 'Denoiser 1: ', est1[:display], '(length ', len(est1), ' error ', error(est1, x)/(n+0.0), ')'
    #print 'Denoiser 2: ', est2[:display], '(length ', len(est2), ' error ', error(est2, x)/(n+0.0), ')'
    print '\n'*5

alphas = [0.01, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99]
epsilons = [0.01, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99]
for i in range(len(alphas)):
    for j in range(len(epsilons)):
        markovSourceDenoise(alphas[i], epsilons[j])            
