'''
Randomization Experiment with DSEG Algorithm for CNVs
(n) Seth Temple (e) sdtemple at uw dot edu
'''

from dsegs import *

i = open('../H3/chm13.chr16.txt', 'r')

reads = []
for line in i:
    line = line.strip()
    eol = len(line) - 1
    t1 = line.find('\t')
    t2 = line.find('\t', t1 + 1, eol)
    reads.append(int(line[(t2 + 1):]))

v = {0:-0.1282, 1:0.5649, 2:1.2581, 3:1.9842}
D = -4
S = [4, 8, 12, 16, 20]

sims = []
for j in range(10):
    r = rdm_readcts(reads, 3, len(reads))
    sims.append([dseg_stats(dseg(r, v, 3, s, D)) for s in S])

rtable = [dseg_stats(dseg(reads, v, 3, s, D)) for s in S]
          
stable = []
for k in range(len(S)):
    row = []
    row.append(min([sim[k][0] for sim in sims]))
    row.append(max([sim[k][1] for sim in sims]))
    row.append(sum([sim[k][2] for sim in sims]) / 10)
    stable.append(row)
  
i.close()
