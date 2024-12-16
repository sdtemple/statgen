'''
PhyloHMM Viterbi Parse Implementation
(n) Seth Temple (e) sdtemple at uw dot edu
'''

from gsci_utility import *
from twostate_hmm import *
from math import log2

# file names for alignment counts
f0 = 'STATE1_anc_rep_counts.txt'
f1 = 'STATE2_codon1_2_counts.txt'
with open(f0, 'r') as o0:
    c0 = {}
    for line in o0:
        line = line.strip()
        align = line[:3]
        count = int(line[4:])
        c0[align] = count
        

with open(f1, 'r') as o1:
    c1 = {}
    for line in o1:
        line = line.strip()
        align = line[:3]
        count = int(line[4:])
        c1[align] = count

# emissions
e0 = log_emissions(align_freqs(c0))
e1 = log_emissions(align_freqs(c1))

# read in alignment observation sequence
alignment = 'ENm010.aln.txt'
a = open(alignment, 'r')
obs = malign_parse(a)
a.close()

# transitions
i = [log2(.95), log2(.05)]
t = log_transitions([(.95, .05),(.10, .90)])

# viterbi algorithm
v = viterbi_parse(obs, i, t, [e0, e1])
hst = hist_states(v[1])
hsg = hist_segments(v[1], obs)
        
