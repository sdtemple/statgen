'''
.fna to .fasta file transformation
(n) Seth Temple (e) sdtemple@uw.edu

Command line: python fna_to_fasta.py input output #lines
'''


from sys import argv
fna = open(str(argv[1]),'r')
fasta = open(str(argv[2]),'w')
fasta.write(fna.readline())
for i in range(int(argv[3])):
    fasta.write(fna.readline().strip('\n'))
fasta.close()
fna.close()
