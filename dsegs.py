'''
DSEGS for Read Count Data
(n) Seth Temple (e) sdtemple at uw dot edu
'''

from random import random

def dseg(reads, v, sep, S, D):
    '''
    Determines D-segments of a sequence
    :param reads: list of read counts
    :param v: scoring dict
    :param sep: censoring int
    :param S: num score threshold
    :param D: num drop threshold
    :return: D-segments
    :rtype: list
    '''
    segs = []
    total = 0
    high = 0
    start = 0 # python indices begin with 0
    end = 0
    n = len(reads)
    for j in range(n):

        # add score
        if(reads[j] >= sep):
            total += v[sep]
        else:
            total += v[reads[j]]

        # possible update    
        if(total >= high):
            high = total
            end = j

        # determine D-segment
        if((total <= 0) or (total <= (high + D)) or (j == (n - 1))):
            if(high >= S):
                segs.append((start, end, high))
            high = 0
            total = 0
            start = j + 1
            end = j + 1

    return(segs)

def dseg_stats(segs):
    '''
    Provide some summary of the D-segments
    :param segs: D-segments
    :return: # D-segments, and min/max scores
    :rtype: tuple
    '''
    if(len(segs) == 0):
        return((-1.00, -1.00, 0.00))
    scores = [item[2] for item in segs]
    return((min(scores), max(scores), len(segs)))

def rdm_readcts(reads, sep, vlen):
    '''
    Generate randomized read count data
    :param reads: list of read counts
    :param sep: censoring int
    :param vlen: length of return vector
    :return: randomized read counts
    :rtype: list
    '''
    n = len(reads)
    c = []
    for i in range(sep):
        c.append(sum([1 for ct in reads if ct == i]))
    c.append(n - sum(c))
    r = []
    for i in range(vlen):
        x = random()
        ct = c[0]
        j = 0
        while(x >= ct / n):
            j += 1
            ct += c[j]
        r.append(j)
    return(r)
            
            








        
