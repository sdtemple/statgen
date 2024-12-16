'''
Genomic Science Utility Functions
(n) Seth Temple (e) sdtemple at uw dot edu
'''

def reverse_gnm(gnm):
    rgnm = gnm[::-1]

    # AT
    rgnm = rgnm.replace('a','A')
    rgnm = rgnm.replace('t','a')
    rgnm = rgnm.replace('A','t')

    # CG
    rgnm = rgnm.replace('c','C')
    rgnm = rgnm.replace('g','c')
    rgnm = rgnm.replace('C','g')
    
    return(rgnm)

def complement_gnm(gnm):
    cgnm = gnm

        # AT
    cgnm = cgnm.replace('a','A')
    cgnm = cgnm.replace('t','a')
    cgnm = cgnm.replace('A','t')

    # CG
    cgnm = cgnm.replace('c','C')
    cgnm = cgnm.replace('g','c')
    cgnm = cgnm.replace('C','g')

    return(cgnm)

def nucl_cts(gnm, ctr = {'a':0, 'c':0, 'g':0, 't':0, 'n':0}):
    for ch in gnm:
        ctr[ch] += 1
    return(ctr)

def nucl_freqs(ctr, frq = {'a':.25, 'c':.25, 'g':.25, 't':.25}):
    '''
    Calculate nucleotide frequencies for a genome
    :param ctr: count dict for 1 strand
    :param frq: frequency dict
    :return: frequency dict
    :rtype: dict
    '''
    total = 0
    for na in ['a', 'c', 'g', 't']:
        total += ctr[na]
        
    frq['a'] = float(ctr['a'] + ctr['t']) / total / 2
    frq['t'] = frq['a']
    frq['c'] = float(ctr['c'] + ctr['g']) / total / 2
    frq['g'] = frq['c']
        
    return(frq)

def align_freqs(cts):
    '''
    Calculate alignment frequencies
    :param cts: dict for alignment counts
    :return: alignment frequencies
    :rtype: dict
    '''
    counts = [val for (key, val) in cts.items()]
    total = sum(counts)
    freqs = {key: float(val) / total for (key, val) in cts.items()}
    freqs = {key: freqs[key] for key in sorted(freqs.keys())}
    return freqs

def malign_parse(file):
    '''
    Set up PhlyoHMM observation sequence
    :param f: file pointer
    :return: alignment observations indexed by position
    :rtype: list
    '''
    obs = []
    while(True):
        header = file.readline().strip()
        org1 = file.readline().strip()
        org2 = file.readline().strip()
        org3 = file.readline().strip()
        file.readline()

        if header == '':
            break

        start = int(header[(header.find(':') + 1):header.find('-')])
        end = int(header[(header.find('-') + 1):])
        org1 = org1[(org1.find('\t') + 1):]
        org2 = org2[(org2.find('\t') + 1):]
        org3 = org3[(org3.find('\t') + 1):]

        for i in range(0, end - start + 1):
            obs.append((start + i, org1[i] + org2[i] + org3[i]))

    return obs
        
        


            
