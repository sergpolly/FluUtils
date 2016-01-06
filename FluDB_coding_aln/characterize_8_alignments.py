import re
import os
import sys
from Bio import Seq
from Bio import SeqIO
from Bio import SeqRecord
# from Bio import Align
from Bio import AlignIO
# import itertools
import numpy as np
import pandas as pd


# potential interesting columns = ['name','seg','pos']
seg_num = 8
# name = 'test1'

# typical fname for segment alignemnt ...
aligned_seg_fname = lambda number: "seg%d.afa"%number

# load coding part alignments for all of the segments ...
segment_aln = {}
for seg_idx in range(1,seg_num+1):
    fname = aligned_seg_fname(seg_idx)
    segment_aln['seg%d'%seg_idx] = AlignIO.read(fname,'fasta')
    # brief alignemnt info
    aln_size = len(segment_aln['seg%d'%seg_idx])
    aln_length = segment_aln['seg%d'%seg_idx].get_alignment_length()
    # print
    # print "segment %d: %d sequences of length %d are aligned (length%%3=%d)"%(seg_idx, aln_size, aln_length, aln_length%3)


# as all of the alignments are loaded ...
# come up with some criteria for loci of interest ...
consensus_genome = []
genome_name_suffix = '_pH1N1'
# let's get some brief statistical info at first ...
for seg in sorted(segment_aln):
    aln = np.array([list(rec) for rec in segment_aln[seg]], np.character)
    aln = pd.DataFrame(aln)
    # replace gaps with None, so pd.isnull is the way to check gaps ...
    aln = aln.where(aln!='-',None)
    # features of interest ...
    descr = aln.describe()
    # count unique top freq
    # descr.loc['top'] - is a consensus sequence
    descr.loc['freq_ratio'] = descr.loc['freq']/descr.loc['count'] # freq_ratio ignores gaps ...
    descr.loc['gaps'] = aln.isnull().sum()
    #
    #
    print
    print "************************** %s **************************"%seg
    print "%d sequences of length %d are aligned (length%%3=%d)"%(aln.shape[0],aln.shape[1],aln.shape[1]%3)
    # print "Brief stat. info on segment %s"%seg
    # number of variable positions ...
    pos_with_var = (descr.loc['freq']<descr.loc['count']).sum()
    print "There are %d position with some variability, overall %.1f%%"%(pos_with_var,pos_with_var*100.0/aln.shape[1])
    # positions with lot of variation ...
    freq_threshold = 0.92
    print "There are %d positions of signif. variability (ref<%.2f)"%((descr.loc['freq_ratio']<freq_threshold).sum(),
                                                                        freq_threshold )
    #
    # the most variable position ...
    var_pos = descr.loc['freq_ratio'].argmin()
    most_var_pos_content = aln[var_pos].value_counts().to_dict()
    print "Most variable position is %d (zero-indexing), val counts:"%var_pos, str(most_var_pos_content).strip('{}').replace('\'','').replace(': ',':')
    #
    # number of gaps and their location ...
    n_gaps = descr.loc['gaps'].sum()
    gaps_pos = descr.columns[descr.loc['gaps']>0].tolist()
    print "There are %d gaps and their positions are "%n_gaps, str(gaps_pos).strip('[]').replace(' ','')
    print "************************************************************"
    #
    #
    ############## consensus sequence output ####################
    #
    seg_consensus = Seq.Seq(''.join(descr.loc['top']))
    seg_consensus = SeqRecord.SeqRecord(seg_consensus,id="%s%s"%(seg,genome_name_suffix),description='')
    consensus_genome.append(seg_consensus)



# output consensus genome here ...
consens_fname = 'pH1N1.fa'
print
print "Writing consensus genome to %s"%consens_fname
SeqIO.write(consensus_genome,consens_fname,'fasta')



































