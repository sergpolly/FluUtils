import re
import os
import sys
from Bio import Seq
from Bio import SeqIO
import pandas as pd
import itertools
import numpy as np
import random
import subprocess as sub

def get_least_gaps_seq(seq_dict,length,side='left'):
	middle = length/2
	min_gaps = length
	min_key  = ''
	# for all sequences check number of gaps in either half and return its id ...
	for seq_id in seq_dict:
		seq_half = seq_dict[seq_id].seq[:middle] if side=='left' else seq_dict[seq_id].seq[middle:]
		num_gaps = seq_half.count('-')
		# reassign the min gaps counts and id in a procedural fashion ...
		if num_gaps < min_gaps:
			min_gaps = num_gaps
			min_key = seq_id
	# return ...
	return (min_key, min_gaps)

# command to clust sequences and get a draft alignment ...
# usearch -cluster_fast seg1.fasta -id 0.993 -centroids nr.fasta -uc clust.uc

path = "/home/venevs/fludb_pH1N1"


if len(sys.argv) < 3:
	print "Call signature is: \"%s msa_fname subs_size\""
msa_fname = sys.argv[1]
subs_size = int(sys.argv[2])


msa = SeqIO.parse(msa_fname,"fasta")
msa = SeqIO.to_dict(msa)

# chosen msa keys ...
chosen_keys = random.sample(msa,subs_size)

# add sequences with the longest UTRs as well ...
alignment_len = len(msa[chosen_keys[0]].seq)

# find sequence with the least gaps in the left half of the sequence ...
# supposedly - longest left-UTR
left_utr_key,_ = get_least_gaps_seq(msa,alignment_len,side='left')
# find sequence with the least gaps in the right half of the sequence ...
# supposedly - longest right-UTR
right_utr_key,_ = get_least_gaps_seq(msa,alignment_len,side='right')


# include those 2 if they are yet in the subsampled alignement ..
if left_utr_key not in chosen_keys:
	chosen_keys += [left_utr_key, ]
if right_utr_key not in chosen_keys:
	chosen_keys += [right_utr_key, ]

# now extract aligned sequences ...
alignment_out = [msa[sid] for sid in chosen_keys]
# output the alignment now ...
tmp_afa_fname = "tmp.afa"
SeqIO.write(alignment_out,tmp_afa_fname,"fasta")
# htm out fname :
out_htm = os.path.basename(msa_fname)+'.htm'
cmd = "mview -in fasta -ruler on -moltype dna -coloring consensus -threshold 60 -consensus on -con_threshold 60 -html head %s > %s"%(tmp_afa_fname,out_htm)
print
print cmd
print
retcode = sub.call(cmd,shell=True)
if retcode == 0:
	print "Complete ..."
else:
	print "mview retcode was %s"%str(retcode)
# # 
# # remove temporary file here ...
# os.remove(tmp_afa_fname)
# print "tmp file removed ..."
# # now make an html alignment using mview ...










































































































