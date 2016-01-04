import re
import os
import sys
from Bio import Seq
from Bio import SeqIO
from Bio import SeqRecord
# import pandas as pd
# import itertools
import numpy as np
# import time


absolute_threshold = 22
delta_coefficient = 0.333
defect_seq_limit = 5


# given an array of integer values, it creates a histogramm (dictionary) ...
def int_hist(array):
    values_dict = dict((val,0) for val in set(array))
    for val in array:
        values_dict[val] += 1
    return values_dict

# hamming distance betwen 2 sequences ...
def hamming(seq1,seq2):
    # sequence length must match ...
    assert len(seq1)==len(seq2)
    return sum(letter1!=letter2 for letter1,letter2 in zip(seq1,seq2))


def get_aligned_pos(sequence,primer,direction='forward'):
    primer_len = len(primer)
    seq_len = len(sequence)
    delta_threshold = primer_len*delta_coefficient
    # offset
    offset = ''.join('_' for char in primer)
    # add offest in the from or in the back, depending on the scan type ...
    seq_check = (offset+sequence) if direction=='forward' else (sequence+offset)
    # position_aligned
    # the range to scan ...
    scan_range = range(0+1, seq_len-1 ) if direction=='forward' else range(seq_len-1,0,-1)
    for frame in scan_range:
        # slide along the reference sequence to find where primer matches the best ...
        mis_left    = hamming( primer, seq_check[frame-1:frame+len(primer)-1] )
        mis_center  = hamming( primer, seq_check[frame+0:frame+len(primer)+0] )
        mis_right   = hamming( primer, seq_check[frame+1:frame+len(primer)+1] )
        # we're looking for a deep 'dip' in the mismatch number ...
        delta_left = mis_center - mis_left
        delta_right = mis_right - mis_center
        # criteria for the perfect alignemnt ...
        crit1 = delta_left  < -delta_threshold
        crit2 = delta_right >  delta_threshold
        crit3 = mis_center  <  absolute_threshold
        #
        if crit1 and crit2 and crit3:
            # print "Position %d does qualify!"%frame
            return (frame-primer_len) if direction=='forward' else frame
    ####################################
    # if nothing happened, then it's bad ...
    print sequence.id
    print str(sequence.seq)
    print direction,primer
    return None


##########################################
def test_mismatch(sequence,primer,direction='forward'):
    mismatch = []
    #
    primer_len = len(primer)
    seq_len = len(sequence)
    delta_threshold = primer_len*delta_coefficient
    # offset
    offset = ''.join('_' for char in primer)
    # add offest in the from or in the back, depending on the scan type ...
    seq_check = (offset+sequence) if direction=='forward' else (sequence+offset)
    # position_aligned
    # the range to scan ...
    scan_range = range(0+1, seq_len-1 ) if direction=='forward' else range(seq_len-1,0,-1)
    for frame in scan_range:
        # slide along the reference sequence to find where primer matches the best ...
        mismatch.append(hamming( primer, seq_check[frame+0:frame+len(primer)+0] ))
    return mismatch

###############################################
def test_mismatch_plot(mmm):
    import matplotlib.pyplot as plt
    plt.plot(mmm,'bo-')
    plt.title( str(sorted(mmm)[:5]) )
    plt.show()


# # given 2 primers - extract only coding part from the segment ...
# primer3 = "ATGGAGAGAATAAAAGAACTGAGAGATCTAATGTCGCAGTCCCGCACTCGCGAGATACTCACTAAGACCACTGTG"
# primer5 = "AAACGAAAACGGGACTCTAGCATACTTACTGACAGCCAGACAGCGACCAAAAGAATTCGGATGGCCATCAATTAG"
def get_primers(fname):
    # simply read 2 lines from the file: 2 lines - 2 primers ...
    with open(fname,'r') as fp:
        left_primer = fp.readline()
        right_primer = fp.readline()
    # unpack primers and stuff ...
    left_name, left_seq  =  left_primer.strip().split(':')
    right_name,right_seq = right_primer.strip().split(':')
    # primer file creator must be aware of the order in primers left first, right next ...
    assert left_name  == 'primer3'
    assert right_name == 'primer5'
    # once all set, just return the primers in the right order ...
    return (left_seq, right_seq)




path = "/home/venevs/fludb_pH1N1"
#
if len(sys.argv)<4:
    print "Call signature: \"%s  in_seq_fname  out_aligned_fname  in_primer_fname\""%sys.argv[0]
    sys.exit(1)
##########################################
seq_fname       = sys.argv[1]
aligned_fname   = sys.argv[2]
primer_fname    = sys.argv[3]
#######################
array_of_seq = SeqIO.parse(seq_fname,"fasta")
array_of_seq = SeqIO.to_dict(array_of_seq)
sorted_sequences_keys = sorted(array_of_seq)
#######################
primer3, primer5 = get_primers(primer_fname)

###############################################################
# CALCULATING START & STOP OF THE CODING REGION ...           #
# ALSO DO SOME OFFSET CALCULATIONS FOR OFFSETS ASSIGNMENT ... #
###############################################################
defect_seq_keys = []
left_off  = []
right_off = []
codlen = []
seqlen = []
for i in sorted_sequences_keys:
    ref_seq = array_of_seq[i]
    ref_seq_len = len(ref_seq)
    start   = get_aligned_pos(ref_seq,primer3,'forward')
    end     = get_aligned_pos(ref_seq,primer5,'backward')
    ####################################
    if (start is None) or (end is None):
        print "Sequence %s is defetive regardless of any futher analysis ..."%i
        defect_seq_keys.append(i)
        # print end-start,start,end
        left_off.append(0)     
        right_off.append(0)    
        codlen.append(0)     
        seqlen.append(0)
    else:     
        # print end-start,start,end
        left_off.append(start-0 +1) # +1 whenever there is delta ...
        right_off.append(ref_seq_len-end +1) # +1 whenever there is delta ...
        codlen.append(end-start +1) # +1 whenever there is delta ...
        seqlen.append(ref_seq_len)
############################
max_left_off  = max(left_off)
max_right_off = max(right_off)


###############################################################################
### DEALING WITH "BAD" SEQUENCES ...
###############################################################################
# assuming all intermediate parts are of equal size ...
the_len_hist = int_hist(codlen)
if len(the_len_hist) > 1:
    print "Length variation of the coding part detected ..."
# simply set the limit on the sequences that can be thrown away,
# and throw them away ...
# sort histogram by values, and get the least populated lengths (the outliers or defects) ...
lens_of_defects = sorted(the_len_hist, key=the_len_hist.get)[:-1]
# most populated length ...
popular_lens = sorted(the_len_hist, key=the_len_hist.get)[-1]
##########################
defect_seq_num = sum(the_len_hist[val] for val in lens_of_defects)
if defect_seq_num < defect_seq_limit:
    print "There are only %d defect sequences. Neglect them and proceed ..."%defect_seq_num
    # find all defective sequecnes ...
    for idx,l in enumerate(codlen):
        if l in lens_of_defects:
            defect_seq_keys.append(sorted_sequences_keys[idx])
else:
    print "There are %d defect sequences. Too much to proceed, but we'll try ..."%defect_seq_num
    # find all defective sequecnes ...
    for idx,l in enumerate(codlen):
        if l in lens_of_defects:
            print l,popular_lens,sorted_sequences_keys[idx]
            defect_seq_keys.append(sorted_sequences_keys[idx])
#############################

# DO ALIGNMENT AT THIS POINT ...
#############################
seq_aligned = []
for i,cid in enumerate(sorted_sequences_keys):
    # avoid defective sequences altogether ...
    if cid not in defect_seq_keys:
        seq = array_of_seq[cid]
        the_left_off = left_off[i]
        the_right_off = right_off[i]
        seq_len = len(seq)
        left_add = ''.join('-' for j in range(max_left_off - the_left_off))
        right_add = ''.join('-' for j in range(max_right_off - the_right_off))
        ################################
        sequence_aligned = Seq.Seq(left_add+str(seq.seq)+right_add)
        seq_aligned.append( SeqRecord.SeqRecord(sequence_aligned,id=seq.id) )
#
# OUTPUT ...
#
SeqIO.write(seq_aligned,aligned_fname,"fasta")



# # plt.plot(mismatch,'bo-')
# # plt.title( str(sorted(mismatch)[:5]) )
# # plt.show()
# # mview -in fasta -ruler on -moltype dna -coloring consensus -threshold 90 -consensus on -con_threshold 90 -html head seqs.afa >yyy.htm

































































































