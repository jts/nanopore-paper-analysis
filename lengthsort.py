from Bio import SeqIO
import sys

def len_and_read_name(a, b):
	length_test = cmp(len(a), len(b))
	if length_test != 0: return length_test
	return cmp(a.id, b.id)

seqs = list(SeqIO.parse(sys.stdin, "fasta"))
seqs.sort(cmp=len_and_read_name, reverse=True)
SeqIO.write(seqs, sys.stdout, "fasta")

