import random
import numpy as np
import math
from pybrain.optimization import CMAES

global kmers
global orig
global noisy
global k
global n
global alphabet
global delta
global rho
global max_del

def allklength(k, alphabet):
    global kmers
    kmers = []
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

def denoiseSequence2(noisy, k, alphabet, deletion_rate, max_del=1, weights=[1]):
    adjust = 1
    if len(weights) != max_del:
        #print('Error: dimension of weights vector is not equal to number of separation lengths.')
        #print('Using default weights.')
        weights = [1]*max_del
    context_hists = []
    for j in range(max_del+1):
        context_del_hist = {}
        for i in range(k+j, len(noisy)-k):
            context_del = noisy[i-k-j:i-j]+noisy[i:i+k]
            deleted = noisy[i-j:i]
            if context_del in context_del_hist:
                context_del_hist[context_del].append(deleted)
            else:
                context_del_hist[context_del] = [deleted]
        context_hists.append(context_del_hist)
    for j in range(1, max_del+1):
        allklength(j, alphabet)
        for i in range(k, len(noisy)-k):
            context = noisy[i-k:i+k]
            context_hist = context_hists[j]
            context_del_hist = context_hists[0]
            if context in context_hist and context in context_del_hist:
                if adjust*(deletion_rate**j)*len(context_hist[context])*weights[j-1] >= (1.0/adjust)*(1-deletion_rate)*len(context_del_hist[context]):
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
        

def denoiseSequence3(noisy, k, alphabet, rho, l=-1, weights=[1]):
    adjust = 1
    if l == -1:
        l = k
    if len(weights) != l:
        print('Error: dimension of weights vector is not equal to number of votes per context.')
        print('Using default weights.')
        weights = [1]*l
    context_left = {} 
    context_right = {}
    for i in range(len(noisy)):
        if i >= k:
            lcontext = noisy[i-k:i]
            if lcontext in context_left:
                context_left[lcontext] += noisy[i]
            else:
                context_left[lcontext] = [noisy[i]]
        if i < len(noisy)-k:
            rcontext = noisy[i+1:i+k+1]
            if rcontext in context_right:
                context_right[rcontext] += noisy[i]
            else:
                context_right[rcontext] = [noisy[i]]
    for i in range(k, len(noisy)-k-1):
        context1 = noisy[i-k+1:i+k+2]
        for a in alphabet:
            vote = 0
            context2 = noisy[i-k+1:i+1]+a+noisy[i+1:i+k+1]
            if a+rcontext[:-1] in context_right and lcontext[1:]+a in context_left:
                for x in range(l):
                    if context1[i-x+1:i-x+k+1] in context_right and context2[i-x+1:i-x+k] in context_right:
                        p1 = sum([int(y == noisy[i-x]) for y in context_right[context1[i-x+1:i-x+k+1]]])
                        p2 = sum([int(y == noisy[i-x]) for y in context_right[context2[i-x+1:i-x+k+1]]])
                        if (1-rho)*p1 <= rho*p2:
                            vote += weights[x]
                for x in range(l):
                    if context1[i-x+1:i-x+k+1] in context_left and context2[i+x-k:i+x] in context_left:
                        p1 = sum([int(y == noisy[i+x]) for y in context_left[context1[i+x-k:i+x]]]) 
                        p2 = sum([int(y == noisy[i+x]) for y in context_left[context2[i+x-k:i+x]]])
                        if (1-rho)*p1 <= rho*p2:
                            vote += weights[x]
                if vote >= l:
                    noisy = noisy[:i+1]+a+noisy[i+1:]
    return noisy

def textDenoise(filename):
    n = 100000
    ks = open(filename, 'r').read()[:n]
    k = int(0.2*math.log(n, 2))
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

def weightWrapper1(weights):
    est1 = denoiseSequence2(noisy, k, alphabet, rho, max_del, weights)
    est = optimalDenoise(noisy, k, alphabet, rho, delta)
    err = error(est1, orig)/(n+0.0)
    return err
    
def weightWrapper2(weights):
    est2 = denoiseSequence3(noisy, k, alphabet, rho, max_del, weights)
    est = optimalDenoise(noisy, k, alphabet, rho, delta)
    err = error(est2, orig)/(n+0.0)
    return err

def generateSequence(n, a):
    p = random.random()
    global k
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
    return x

def markovSourceDenoise(transition_prob, deletion_rate, weighted=False):
    global delta
    delta = transition_prob
    global rho
    rho = deletion_rate
    global max_del
    max_del = 2
    global n
    n = 10000
    display = 50
    global alphabet
    alphabet = ['+', '-']
    global orig
    orig = generateSequence(n, delta)
    global noisy
    noisy = deletionChannel(orig, rho)
    est1 = denoiseSequence2(noisy, k, alphabet, rho, max_del)
    #est2 = denoiseSequence3(noisy, k, alphabet, rho)
    est = optimalDenoise(noisy, k, alphabet, rho, delta)
    print 'Setting: delta = ', delta, ', rho = ', rho 
    print 'Original: ', orig[:display], '(length ', len(orig), ' error ', error(orig, orig)/(n+0.0), levenshtein(orig, orig)/(n+0.0), ')'
    print 'Noisy: ', noisy[:display], '(length ', len(noisy), ' error ', error(noisy, orig)/(n+0.0), levenshtein(noisy, orig)/(n+0.0), ')'
    print 'Denoised: ', est[:display], '(length ', len(est), ' error ', error(est, orig)/(n+0.0), levenshtein(est, orig)/(n+0.0), ' )'
    print 'Denoiser 1: ', est1[:display], '(length ', len(est1), ' error ', error(est1, orig)/(n+0.0), levenshtein(est1, orig)/(n+0.0), ')'
    #print 'Denoiser 2: ', est2[:display], '(length ', len(est2), ' error ', error(est2, orig)/(n+0.0), ')'
    if weighted:
        l1 = CMAES(weightWrapper1, [1])
        l1.minimize = True
        l1.maxEvaluations = 200
        opt_vals1 = l1.learn()
        weights1 = opt_vals1[0]
        err1 = opt_vals1[1]
        #l2 = CMAES(weightWrapper2, [1]*k)
        #l2.minimize = True
        #l2.maxEvaluations = 200
        #opt_vals2 = l2.learn()
        #weights2 = opt_vals2[0]
        #err2 = opt_vals2[1]
        print(weights1)
        weighted_est1 = denoiseSequence2(noisy, k, alphabet, rho, max_del, weights1)
        #weighted_est2 = denoiseSequence3(noisy, k, alphabet, rho, k, weights2)
        print 'Weighted Denoiser 1: ', weighted_est1[:display], '(length ', len(weighted_est1), ' error ', err1, ')'
        #print 'Weighted Denoiser 2: ', weighted_est2[:display], '(length ', len(weighted_est2), ' error ', err2, ')'
    print '\n'*5

def main():
    alphas = [0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99]
    epsilons = [0.01, 0.1, 0.2, 0.3, 0.4]
    for a in alphas:
        for e in epsilons:
            markovSourceDenoise(a, e, True)

main()
