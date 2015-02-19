import sys
from Bio import SeqIO

dna_alphabet = ['A', 'C', 'G', 'T']

# Generate a list of all k-mers
def generate_mers(alphabet, k):

    # initialze kmer list with an empty string
    kmers = list()
    kmers.append("")

    for i in range(0, k):

        # for every string in the output list,
        # generate a new string one base longer for
        # every symbol in the alphabet
        new_kmers = list()
        for o in kmers:
            for j in range(0, len(alphabet)):
                new_kmers.append(o + alphabet[j])
        kmers = new_kmers

    kmers.sort()
    return kmers

# Generate a map from k-mer -> lexicographic rank
def generate_mers_rank_map(alphabet, k):

    kmers = generate_mers(alphabet, k)

    # Build dictionary
    rank_dict = dict()
    for (i, o) in enumerate(kmers):
        rank_dict[o] = i
    return rank_dict

# Generate all of the kmers of the string
def str2kmers(s, k):
    out = []
    for i in xrange(0, len(s) - k + 1):
        out.append(s[i:i+k])
    return out

# reverse complement a sequence
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

def revcomp(seq):
    reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
    return reverse_complement

#
# Main
#
K = 5
kmer_set = generate_mers(dna_alphabet, K)
rank_map = generate_mers_rank_map(dna_alphabet, K)

counts_by_name = dict()

for fn in sys.argv[1:]:
    for rec in SeqIO.parse(open(fn), "fasta"):
        fwd_name = rec.name + ".fwd"
        rc_name = rec.name + ".rc"
        
        counts_by_name[fwd_name] = [0] * (4 ** K)
        counts_by_name[rc_name] = [0] * (4 ** K)
        
        for ki in xrange(0, len(rec.seq) - K + 1):
            kmer = rec.seq[ki:ki+K]
            
            if 'n' in kmer or 'N' in kmer:
                continue

            fwd_rank = rank_map[str(kmer)]
            rc_rank = rank_map[revcomp(str(kmer))]

            counts_by_name[fwd_name][fwd_rank] += 1
            counts_by_name[rc_name][rc_rank] += 1

datasets = sorted(counts_by_name.keys())

# header
print "\t".join(["kmer"] + datasets)

# counts per kmer
for kmer in kmer_set:
    
    out = list()
    out.append(kmer)
    rank = rank_map[kmer]  

    for d in sorted(datasets):
        out.append(str(counts_by_name[d][rank]))
    print "\t".join(out)
