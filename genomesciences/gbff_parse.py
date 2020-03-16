'''
.gbff Parse Utility
(n) Seth Temple (e) sdtemple at uw dot edu
'''

# global
fstrand = []
cstrand = []

def to_features(file):
    '''
    Move file pointer to features list
    :param file: file pointer
    :rtype: None 
    '''
    line = file.readline()
    while(line[:8] != 'FEATURES'):
        line = file.readline()
    return(None)

def is_feature(line):
    '''
    Identify if line is the start of a feature
    :param line: str
    :rtype: bool
    '''
    bl = line[5] != ' '
    return(bl)

def to_feature(file):
    '''
    Move file pointer to a feature
    :param file: file pointer
    :return: start line for a feature
    :rtype: str
    '''
    line = file.readline()
    while not is_feature(line):
        line = file.readline()
        if line[:6] == 'ORIGIN':
           break
    return(line)

def is_cds(line):
    '''
    Identify if a feature is a coding sequence
    :param line: str
    :rtype: bool
    '''
    bl = line[5:8] == 'CDS'
    return(bl)

def to_cds(file):
    '''
    Move file pointer to a CDS feature
    :param file: file pointer
    :return: start line for a CDS feature
    :rtype: str
    '''
    line = file.readline()
    while not is_cds(line):
          line = to_feature(file)
          if line[:6] == 'ORIGIN':
              break
    return(line)

def get_cds(file, line):
    '''
    Get the coding sequence
    :param file: file pointer
    :param line: start line for a CDS feature
    :return: entire coding sequence
    :rtype: str
    '''
    seq = ''
    line = (line.strip()[3:]).strip()
    while(line[0] != '/'):
        seq += line
        line = file.readline().strip()
    return(seq)
    
def assign_strand(seq):
    '''
    Assign coding sequence to a strand
    :param seq: str coding sequence
    :rtype: None
    '''

    # disregard ambiguous starts and ends
    if '<' in seq:
        return(None)
    if '>' in seq:
        return(None)

    # cds on forward strand
    if '(' not in seq:
        period = seq.find('.')
        start = int(seq[:period]) - 1
        end = int(seq[(period + 2):]) - 1
        fstrand.append(write_splice(start, end, seq))

    # long cds on reverse strand
    elif 'complement(join(' in seq:
        seq = seq[16:-2]
        period = seq[::-1].find('.')
        comma = seq[::-1].find(',')
        if period > comma:
            start = int(seq[-comma:]) - 1
            end = start
        else:
            # only subseqs of size 1 in sequence
            if period == -1:
                end = int(seq[-comma:]) - 1
                start = end
            # most common
            else:
                end = int(seq[-period:]) - 1
                start = int(seq[-comma:-(period + 2)]) - 1
        cstrand.append(write_splice(start, end, seq, -1))

     # cds on reverse strand       
    elif 'complement(' in seq:
        seq = seq[11:-1]
        period = seq[::-1].find('.')
        end = int(seq[-period:]) - 1
        start = int(seq[:-(period + 2)]) - 1
        cstrand.append(write_splice(start, end, seq, -1))

    # long cds on forward strand     
    elif 'join(' in seq:
        seq = seq[5:-1]
        period = seq.find('.')
        comma = seq.find(',')
        if period > comma:
            start = int(seq[:comma]) - 1
            end = start
        else:
            # only subseqs of size 1 in sequence
            if period == -1:
                start = int(seq[:comma]) - 1
                end = start
            # most common
            else:
                start = int(seq[:period]) - 1
                end = int(seq[(period + 2):comma]) - 1
        fstrand.append(write_splice(start, end, seq))

    return(None)
        

def write_exon(i, j, ct = 0, d = 1):
    '''
    Subroutine to write partial exon exon
    :param i: int
    :param j: int
    :param ct: int counter
    :return: exon
    :rtype: list
    '''
    l = []
    if ct > 0:
        if d == -1:
            for k in range(j, i - 1, d):
                if ct == 10:
                    break
                l.append(k)
                ct += 1
        else:
            for k in range(i, j + 1, d):
                if ct == 10:
                    break
                l.append(k)
                ct += 1
    
    elif ct == 0:
        if d == -1:
            l.append(j)
            for k in range(j - 1, i - 1, d):
                if ct == 10:
                    break
                l.append(k)
                ct += 1
        else:
            l.append(i)
            for k in range(i + 1, j + 1, d):
                if ct == 10:
                    break
                l.append(k)
                ct += 1

    elif ct < 0:
        if d == -1:
            ct = 0
            for k in range(j, i - 1, d):
                if ct == 10:
                    break
                l.append(k)
                ct += 1
        else:
            ct = 0
            for k in range(i, j + 1, d):
                if ct == 10:
                    break
                l.append(k)
                ct += 1
    
            
    return(l)

def write_exons(i, j, seq, d = 1):
    '''
    Routine to write partial exonic sequence
    :param i: int
    :param j: int
    :param seq: str coding sequence
    :return: partial exonic sequence
    :rtype: list
    '''
    ct = (j - i) if (j - i) > 0 else -1
    asplice = write_exon(i, j, 0, d)
    while(ct < 10):
        
        if d == -1:
            comma = seq[::-1].find(',')
            seq = seq[:-(comma + 1)]
            period = seq[::-1].find('.')
            comma = seq[::-1].find(',')
            if period > comma and comma != -1:
                end = int(seq[-comma:]) - 1
                start = end
            else:
                # only number in sequence
                if period == -1 and comma == -1:
                    start = int(seq) - 1
                    end = start
                # only subseqs of size 1 in sequence
                elif period == -1 and comma != -1:
                    end = int(seq[-comma:]) - 1
                    start = end
                # most common
                else:
                    end = int(seq[-period:]) - 1
                    if comma == -1:
                        start = int(seq[:-(period + 2)]) - 1
                    else:
                        start = int(seq[-comma:-(period + 2)]) - 1
                
        else:
            comma = seq.find(',')
            seq = seq[(comma + 1):]
            period = seq.find('.')
            comma = seq.find(',')
            if period > comma and comma != -1:
                start = int(seq[:comma]) - 1
                end = start
            else:
                # one number in sequence
                if period == -1 and comma == -1:
                    start = int(seq) - 1
                    end = start
                # only subseqs of size 1 in sequence
                elif period == -1 and comma != -1:
                    start = int(seq[:comma]) - 1
                    end = start
                # most common
                else:
                    start = int(seq[:period]) - 1
                    if comma == -1:
                        end = int(seq[(period + 2):]) - 1
                    else:
                        end = int(seq[(period + 2):comma]) - 1
              
        asplice += write_exon(start, end, ct, d)
        ct += (end - start + 1)
        
    return(asplice)
        

def before_splice(i, j, d = 1):
    '''
    Routine to write 10 positions before splice site
    :param i: int
    :return: sequence before splice site
    :rtype: list
    '''
    l = []
    if d == -1:
        for k in range(1,11):
            l.append(j + k)
    else:
        for k in range(1,11):
            l.append(i - k)
    l.reverse()
    return(l)


def write_splice(i, j, seq, d = 1):
    '''
    Routine to write 21 positions surrounding splice site
    :param i: int
    :param j: int
    :param seq: str coding sequence
    :return: 21 positions surrounding splice site
    :rtype: list
    '''
    bsplice = before_splice(i, j, d)
    asplice = write_exons(i, j, seq, d)
    return(bsplice + asplice)

def to_origin(file):
    '''
    Parse .gbff until ORIGIN line
    :param file: file pointer
    :rtype: None
    '''
    while(True):
        line = to_cds(file)
        if line[:6] == 'ORIGIN':
            break
        else:
            assign_strand(get_cds(file, line))

    return(None)

def get_nucl(line):
    '''
    Parse .gbff lines for partial genome
    :param line: str
    :return: partial genome
    :rtype: str
    '''
    skip = line.find(' ')
    line = line[(skip + 1):]
    line = line.replace(' ','')
    return(line)

def get_gnm(file):
    '''
    Parse .gbff file for genome
    :param file: file pointer
    :return: genome
    :rtype: str
    '''
    gnm = ''
    line = get_nucl(file.readline().strip())
    while line != '//':
        gnm += line
        line = get_nucl(file.readline().strip())
    return(gnm)

def count_sites(fss, rss, fs, rs, n = 21):
    '''
    Count nucleotides at each site
    :param fss: forward strand sites
    :param rss: reverse strand sites
    :param fs: forward strand
    :param rs: reverse strand
    :param n: int
    :return: counts of nucleotides in target
    :rtype: dict
    '''
    ctr = [{'a':0, 'c':0, 'g':0, 't':0, 'n':0} for i in range(n)]

    for site in fss:
        for i in range(len(site)):
            if i > 20:
                print(i)
                print(site)
            ctr[i][fs[site[i]]] += 1

    for site in rss:
        for i in range(len(site)):
            if i > 20:
                print(i)
                print(site)
            ctr[i][rs[site[i]]] += 1

    return(ctr)

def freq_site(sctr):
    '''
    Calculate nucleotide frequencies at a site
    :param sctr: counts dictionary
    :return: frequencies dictionary
    :rtype: dict
    '''
    frq = {'a':0, 'c':0, 'g':0, 't':0}
    total = float(sctr['a'] + sctr['c'] + sctr['g'] + sctr['t'])
    frq['a'] = sctr['a'] / total
    frq['c'] = sctr['c'] / total
    frq['g'] = sctr['g'] / total
    frq['t'] = sctr['t'] / total
    return(frq)

from math import log2
def wt_site(sfrq, gfrq, null = -99):
    '''
    Calculate LLR weights at a site
    :param sfrq: site specific frequencies dictionary
    :param gfrq: genome frequencies dictionary
    :return: weights dictionary
    :rtype: dict
    '''
    wt = {'a':0, 'c':0, 'g':0, 't':0, 'n':0}
    wt['a'] = log2(sfrq['a']) - log2(gfrq['a']) if sfrq['a'] > 0 else null
    wt['c'] = log2(sfrq['c']) - log2(gfrq['c']) if sfrq['c'] > 0 else null
    wt['g'] = log2(sfrq['g']) - log2(gfrq['g']) if sfrq['g'] > 0 else null
    wt['t'] = log2(sfrq['t']) - log2(gfrq['t']) if sfrq['t'] > 0 else null
    return(wt)

def score_site(site, strand, wts):
    '''
    Calculate the score for a possible splice site
    :param site: list of positions
    :param strand: str
    :param wts: list of weights dictionaries
    :return: score of site
    :rtype: float
    '''
    site = [strand[i] for i in site]
    scores = [wts[i][site[i]] for i in range(len(site))]
    return(sum(scores))
    




