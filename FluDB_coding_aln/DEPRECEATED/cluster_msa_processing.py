import re
import os
import sys
from Bio import Seq
from Bio import SeqIO
import pandas as pd
import itertools
import numpy as np

path = "/home/venevs/fludb_pH1N1"
# seg5-GenomicFastaResults.fasta
# seg7.afa

centroids = SeqIO.parse(os.path.join(path,"nr.fasta"),"fasta")
clust_msa = SeqIO.parse(os.path.join(path,"seg1.afa"),"fasta")

# just check how centroids are aligned in the clust_msa:
# select centroids and find them in clust_msa -> extract their alignemnt ...


centroids = SeqIO.to_dict(centroids)
clust_msa = SeqIO.to_dict(clust_msa)


# simple test to check if all ids from centroids are present in clust_msa ...
all_match = all([(sid in clust_msa) for sid in centroids])
if not all_match:
	raise ValueError("All ids from centroids must match ids in clust_msa ...")


# now extract aligned centroid sequences ...
centroids_aligned = [clust_msa[sid] for sid in centroids]

# output the alignment now ...
SeqIO.write(centroids_aligned,"centrs.afa","fasta")



# # now let's gather the statistics about each of the segments ...
# # segments with lots of non-ATGC letters ...
# # countATGC = {}
# # length = {}
# data = {}
# for seg in segments:
#     get_ATGC = lambda seq: seq.count('A')+seq.count('T')+seq.count('G')+seq.count('C')
#     get_strain = lambda seqrec: seqrec.description.strip().split('|')[1]
#     countATGC = []
#     length = []
#     strain = []
#     for seqrec in segments[seg]:
#         countATGC.append(get_ATGC(seqrec.seq))
#         length.append(len(seqrec))
#         strain.append(get_strain(seqrec))
#         #
#     data[seg] = pd.DataFrame({"atgc":countATGC,"len":length,"strain":strain})
#     data[seg]["not_atgc"] = data[seg]["len"] - data[seg]["atgc"]

# filtered_index = {}
# length_threshold = 0.9
# #
# print "filtering length cutoff %.2f"%length_threshold

# for seg in segments:
#     filtered_index[seg] = np.asarray(data[seg][(data[seg]["not_atgc"]==0)&(data[seg]["len"]>length_threshold*data[seg]["len"].median())].index)


# # dat = pd.DataFrame({})
# for seg in sorted(segments):
#     msg = "%s: total_seq=%d filtered_seq=%d disregarded_seq=%d"%(seg,len(segments[seg]),filtered_index[seg].size,len(segments[seg]) - filtered_index[seg].size)
#     print msg


# for seg in segments:
#     out_fname = os.path.join(path,"%s.fasta"%seg)
#     SeqIO.write( (segments[seg][idx] for idx in filtered_index[seg]), out_fname, "fasta" )













































































































