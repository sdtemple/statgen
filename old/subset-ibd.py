'''
Filter .ibd file to diploid subset
Seth Temple, sdtemple@uw.edu
'''

#import gzip

### inputs ###
from sys import argv
ibdfile = argv[1]
max_size = int(argv[2])
outfile = argv[3]

### output ###
out = open(outfile, 'w')

### procedure ###
subset = set()
max_size = 100
set_size = 0
with open(ibdfile, 'r') as f:
#with gzip.open(ibdfile, 'r') as f:
    for line in f:
        ab = line.split('\t')
        a = ab[0]
        b = ab[2]
        if a not in subset:
            if set_size < max_size:
                subset.add(a)
                set_size += 1
        if b not in subset:
            if set_size < max_size:
                subset.add(b)
                set_size += 1
        if (b in subset) and (a in subset):
            out.write(line)
            out.write('\n')
out.close()