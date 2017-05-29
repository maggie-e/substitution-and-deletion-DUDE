import random
import numpy as np
import numpy.linalg as LA
import matplotlib.pyplot as plt
from itertools import chain

kmers = []
read_length = 50
num_reads = 10000

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

def calculateMatrix(input_sequence, k, alphabet):
    matrix_dim = len(alphabet)**k;
    transition_matrix = np.zeros((matrix_dim, matrix_dim))
    allklength(k, alphabet)
    sequence_map = {}
    for i in range(len(kmers)):
        sequence_map[kmers[i]] = i
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

def denoiseSequence2(input_sequence, k, alphabet, deletion_rate):
    noisy = deletionChannel(input_sequence, deletion_rate)
    context_hist = {}
    context_del_hist = {}
    for i in range(k, len(noisy)-k+1):
        if i < len(noisy)-k:
            context = noisy[i-k:i]+noisy[i+1:i+k+1]
            if context in context_hist:
                context_hist[context] += noisy[i]
            else:
                context_hist[context] = [noisy[i]]
        context_del = noisy[i-k:i+k]
        if context_del in context_del_hist:
            context_del_hist[context_del] = context_del_hist[context_del] + 1
        else:
            context_del_hist[context_del] = 1
    for i in range(k, len(noisy)-k):
        context = noisy[i-k:i+k]
        if context in context_hist and context in context_del_hist:
            if deletion_rate*len(context_hist[context]) >= (1-deletion_rate)*context_del_hist[context]:
                ml = alphabet[0]
                p = 0
                for a in alphabet:
                    new_p = sum(context_hist[context] == a)/(len(context_hist[context])+0.0)
                    if new_p > p:
                        p = new_p
                        ml = a
                noisy = noisy[:i]+ml+noisy[i:]
    return input_sequence, noisy

######################################### Deletion and substitution channel ############################################

def errorrate_delsub(original,estimate):
    pass

def count_symbols(sequence,alphabet):
    # count the number of symbols in the sequence and returns a LIST of their frequencies
    sequence = np.array(sequence)
    counts = np.zeros(len(alphabet))
    for i,l in enumerate(alphabet):
        counts[i] = sum([int(x == l) for x in sequence])
    return counts

def denoiseSequence3(noisy,k, alphabet, erasedsymbol, ppi, loss): #'noisy' must be a string
    context_hist = {}
    n = len(noisy)
    for i in range(k, len(noisy)-k+1):
        if i < len(noisy)-k:
            context = noisy[i-k:i]+noisy[i+1:i+k+1]
            if context in context_hist:
                context_hist[context] += noisy[i]
            else:
                context_hist[context] = [noisy[i]]
            context_del = noisy[i-k:i+k]
            if context_del in context_hist:
                context_hist[context_del] += erasedsymbol
            else:
                context_hist[context_del] = [erasedsymbol]
    
    # Create a dictionary that tells you what symbol to insert/replace for each context
    new_context_hist = {}
    newalphabet = [i for i in alphabet]
    newalphabet.extend([erasedsymbol])
    for context, symbolswithcontext in context_hist.iteritems():
        symbollist = list(symbolswithcontext)
        values = count_symbols(symbollist, newalphabet)
        new_context_hist[context] = values
    #print new_context_hist
    newseq = [noisy[:k]]
    # Now, make 'noisy' into a list so it's easier to work with
    noisy = list(noisy)
    noisy_mod = [[i,erasedsymbol] for i in noisy[k:len(noisy)-k]] #Insert 'nothing' symbols at alternate spots, leaving the first k and last k symbols alone
    noisy_mod = chain(noisy[:k],[erasedsymbol],noisy_mod,noisy[len(noisy)-k:])
    noisy_mod = list(chain.from_iterable(noisy_mod))
    for j in range(k,2*n-3*k+1):
        if (j+1-k)%2 == 0: # It's a symbol in the original
            jj = (j+1-k)/2+k-1
            # Consider if symbol jj should be substituted by looking at the context
            context = ''.join(noisy[jj-k:jj]+noisy[jj+1:jj+k+1])
            emp_probs = new_context_hist[context]
            front = np.dot(emp_probs.T,LA.inv(ppi))
            storevals = np.zeros(len(newalphabet))
            for i in range(len(newalphabet)):
                storevals[i] = np.dot(front,np.multiply(loss[:,i],emp_probs))
            newseq.extend(newalphabet[np.argmax(storevals)])
        else: # It's an erased symbol
            jj = (j+2-k)/2+k-1
            context = ''.join(noisy[jj-k:jj+k])
            emp_probs = new_context_hist[context]
            front = np.dot(emp_probs.T,LA.inv(ppi))
            storevals = np.zeros(len(newalphabet))
            for i in range(len(newalphabet)):
                storevals[i] = np.dot(front,np.multiply(loss[:,i],emp_probs))
            newseq.extend(newalphabet[np.argmax(storevals)])
    newseq.extend(noisy[len(noisy)-k:])
    estimate = [i for i in newseq if i!=erasedsymbol]
    estimate = "".join(estimate)
    return estimate
# ################################################## Testing the erasure DUDE ####################################################

# n = 10000
# rho = 0.01

# alphabet = ['a', 'c', 't', 'g']
# input = open('humangenome.fasta').read()[:n]
# #input = takeReads(input)
# # erasure_denoised = erasureDenoise(input, k, alphabet, rho) 
# # print errorrate(input, erasure_denoised)
# # #input, denoised = denoiseSequence2(input, k, alphabet, rho)

################################################# Testing the substitution and deletion DUDE ####################################################
# Generate the input sequence with a Markov process
alpha = 0.1 # Probability of the next symbol being the opposite of the current one
# The 0th entry in W,N,X,a,b are dummies. There are n timesteps and the first depends on X[0]. (So X has n+1 entries in total)
n=100
X=np.zeros(n+1) # Original
Y=np.zeros(n+1) # Measured 
W=np.random.uniform(0,1,n+1)
newW = np.zeros(len(W))
newW[W<=alpha]=1
newW[W>alpha]=-1
W=newW
X[0] = -1
currstate = X[0]
pd = 0.1
ps = 0.1
N=np.random.uniform(0,1,n+1)
N = [0 if i<=pd else -1 if pd<i and i<=pd+ps else 1 if i>pd+ps else i for i in N]
for ii in range(1,n+1):
    X[ii]=currstate*W[ii]
    currstate = X[ii]
    Y[ii]=currstate*N[ii]
# Make the sequence generated by the Markov process into a string
X = X[1:]
Y = Y[1:]
erasedsymbol = 'N'
X = ['a' if i==-1 else 'b' for i in X]
str_X = ''.join(X)
Y = ['a' if i==-1 else 'b' if i==1 else '' for i in Y]
str_Y = ''.join(Y)
# Apply the DUDE
k = 2
ll = 1
lossmat = np.array([[0,1,(1+ll)/2.],[1,0,(1+ll)/2.],[ll,ll,0]])
transitionmat = np.array([[1-ps-pd,ps,pd],[ps,1-ps-pd,pd],[0,0,1]])
alphabet = ['a','b']
Yest = denoiseSequence3(str_Y,k,alphabet,erasedsymbol, transitionmat,lossmat)
print 'Original', str_X
print 'Corrupted', str_Y
print 'Estimate from DUDE', Yest
