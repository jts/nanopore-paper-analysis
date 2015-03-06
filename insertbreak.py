from Bio import SeqIO
import sys

breakpos = int(sys.argv[2])

rec = list(SeqIO.parse(open(sys.argv[1]), "fasta"))[0]

SeqIO.write([rec[0:breakpos]], sys.stdout, "fasta")
contigbreak = rec[breakpos:]
contigbreak.id = "Break"
SeqIO.write([contigbreak], sys.stdout, "fasta")
	
