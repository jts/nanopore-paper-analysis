from Bio import SeqIO
from Bio.SeqUtils import GC
import sys

for rec in SeqIO.parse(open(sys.argv[1]), "fasta"):
        for n in xrange(0, len(rec), 1000):
                print n, GC(rec[n:n+1000].seq)
