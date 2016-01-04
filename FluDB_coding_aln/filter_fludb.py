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


segments = {}

for idx in range(8):
    seg_name = "seg%d"%(idx+1)
    segments[seg_name] = list(SeqIO.parse(os.path.join(path,"%s-GenomicFastaResults.fasta"%seg_name),"fasta"))


# now let's gather the statistics about each of the segments ...

# segments with lots of non-ATGC letters ...
# countATGC = {}
# length = {}
data = {}
for seg in segments:
    get_ATGC = lambda seq: seq.count('A')+seq.count('T')+seq.count('G')+seq.count('C')
    get_strain = lambda seqrec: seqrec.description.strip().split('|')[1]
    countATGC = []
    length = []
    strain = []
    for seqrec in segments[seg]:
        countATGC.append(get_ATGC(seqrec.seq))
        length.append(len(seqrec))
        strain.append(get_strain(seqrec))
        #
    data[seg] = pd.DataFrame({"atgc":countATGC,"len":length,"strain":strain})
    data[seg]["not_atgc"] = data[seg]["len"] - data[seg]["atgc"]

filtered_index = {}
length_threshold = 0.9
#
print "filtering length cutoff %.2f"%length_threshold

for seg in segments:
    filtered_index[seg] = np.asarray(data[seg][(data[seg]["not_atgc"]==0)&(data[seg]["len"]>length_threshold*data[seg]["len"].median())].index)


# dat = pd.DataFrame({})
for seg in sorted(segments):
    msg = "%s: total_seq=%d filtered_seq=%d disregarded_seq=%d"%(seg,len(segments[seg]),filtered_index[seg].size,len(segments[seg]) - filtered_index[seg].size)
    print msg


for seg in segments:
    out_fname = os.path.join(path,"%s.fasta"%seg)
    SeqIO.write( (segments[seg][idx] for idx in filtered_index[seg]), out_fname, "fasta" )













































































































