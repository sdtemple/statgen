'''
.fna to .fasta file transformation
(n) Seth Temple (e) sdtemple@uw.edu

Command line python fasta_to_fna.py input output
'''

from sys import argv
fasta = open(argv[1],'r')
fna = open(argv[2],'w')

# header
fna.write(fasta.readline())

# read and write char by char
i = 0
while True:
    c = fasta.read(1)
    if not c:
        break
    fna.write(c)
    i += 1
    if ((i % 70) == 0):
        fna.write('\n') # 70 chars per line

fna.close()
fasta.close()
