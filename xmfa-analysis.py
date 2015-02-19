import sys
from collections import defaultdict
from Bio import SeqIO

def compare(ref, asm):
    
    n_ins = 0
    n_del = 0
    n_sub = 0
    n_match = 0
    n_col = 0
    
    prev_ref_base = -1
    curr_ref_base = int(ref.name.split(":")[1].split("-")[0])

    for (rbase, abase) in zip(ref.seq, asm.seq):

        # update mutation counts
        is_sub = (abase != rbase and rbase != '-' and abase != '-')
        is_ins = (rbase == '-' and abase != '-')
        is_del = (abase == '-' and rbase != '-')
        is_error = is_sub or is_ins or is_del

        if DEBUG and is_error:
            print 'ERROR', curr_ref_base, rbase, abase
        
        n_ins += is_ins
        n_del += is_del
        n_sub += is_sub
        n_match += not is_sub and not is_ins and not is_del

        n_col += 1

        bin_idx = int(curr_ref_base / bin_size)
        base_count[bin_idx] += curr_ref_base != prev_ref_base
        error_count[bin_idx] += is_error

        prev_ref_base = curr_ref_base

        # update reference base
        curr_ref_base += rbase != '-'
        
    n_total = n_sub + n_ins + n_del
    if DEBUG:    
        print "subs: %d (%.2f) ins: %d (%.2f) del: %d (%.2f) total: %d (%.2f) \n" % (n_sub, n_col / n_sub, n_ins, n_col / n_ins, n_del, n_col / n_del, n_total, n_col / n_total)

    return n_match, n_sub, n_del, n_ins
ref_records = list()
assembly_records = list()

bin_size = 10000
base_count = defaultdict(int)
error_count = defaultdict(int)
DEBUG = False

for (ri, rec) in enumerate(SeqIO.parse(open(sys.argv[1]), "fasta")):
    if "NC" in rec.description:
	ref_records.append(rec)
    else:
        assembly_records.append(rec)

max_records = max(len(ref_records), len(assembly_records))

total_match = 0
total_sub = 0
total_del = 0
total_ins = 0

for ri in range(0, max_records):
    if ri >= len(ref_records) or ri >= len(assembly_records):
        continue
    
    ref_record = ref_records[ri]
    assembly_record = assembly_records[ri]
    if len(ref_record.seq) != len(assembly_record.seq):
        continue
    
    if DEBUG:
        print 'Comparing', ref_record.name, 'to', assembly_record.name
        print len(ref_record.seq), len(assembly_record.seq)

    (m, s, d, i) = compare(ref_record, assembly_record)
    total_match += m
    total_sub += s
    total_del += d
    total_ins += i

accuracy = float(total_match) / (total_match + total_sub + total_del + total_ins)

sys.stderr.write("SUMMARY Match: %d Mismatch: %d Deletion: %d Insertion: %d Accuracy: (%.3f)\n" % (total_match, total_sub, total_del, total_ins, accuracy))

for bi in sorted(base_count.keys()):
        print "\t".join([str(x) for x in [bi * bin_size, (bi + 1) * bin_size, base_count[bi], error_count[bi], float(error_count[bi]) / base_count[bi]]])

