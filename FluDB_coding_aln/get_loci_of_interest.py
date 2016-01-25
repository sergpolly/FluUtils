# import re
import os
import sys
# from Bio import Seq
# from Bio import SeqIO
from Bio import SeqRecord
# from Bio import Align
from Bio import AlignIO
# import itertools
import numpy as np
import pandas as pd
#
# add the location of the flu_module ...
sys.path.append("../RNA_SNP_association")
import flu_module as flu
#
######### GENETIC CODE STUFF ######################
from Bio.Data import CodonTable
genetic_code = CodonTable.standard_dna_table.forward_table
stop_codons = dict([ (codon,'*') for codon in CodonTable.standard_dna_table.stop_codons ])
# make genetic code contain both 61 aacid codons and 3 stop codons ...
genetic_code.update(stop_codons)
###################################################
#
######### KD HYDROPHOBICITY SCALE #################
# Kyte-Doolittle hydrophobicity scale ...
from Bio.SeqUtils import ProtParamData
KD = ProtParamData.kd
KD['*'] = 25.0
###################################################
#
# # columns = ['name','seg','pos','description']
# seg_num = 8
# name = 'test1'
#
########################################################
# reference loading ...
ref_fname = "./pH1N1_coding_dat/pH1N1.fa"
orf_fname = "./pH1N1_coding_dat/pH1N1_noPB1F.orf"
ph1n1 = flu.influenza(ref_fname, orf_fname)
#########################################################

# typical fname for segment alignemnt ...
aligned_seg_fname = lambda number: "seg%d.afa"%number

# load coding part alignments for all of the segments ...
segment_aln = {}
for seg_idx in range(1,ph1n1.segnum+1):
    fname = aligned_seg_fname(seg_idx)
    segment_aln['seg%d'%seg_idx] = AlignIO.read(fname,'fasta')
    # brief alignemnt info
    aln_size = len(segment_aln['seg%d'%seg_idx])
    aln_length = segment_aln['seg%d'%seg_idx].get_alignment_length()
    # print
    # print "segment %d: %d sequences of length %d are aligned (length%%3=%d)"%(seg_idx, aln_size, aln_length, aln_length%3)
#
#######################################
# we would need to import Standard Codon Table from Bio
# we would need the flu module to load ORFs
# some extra function would be needed to understand if a loci/SNP is synonymous/nonsyn/etc.
# we would need something like KD scale just to guess how drastic the amino acid substitution is ...
#######################################
# flu_module ORF system has been updated and is working now.
#######################################
# genetic code and KD loaded ...
#######################################
#  this function is working fine ...
def get_loci_orf_assoc(orf,loci):
    """checks if a given position is in any of the ORFs. Returns
    a list of the hosting ORFs with corresponding coords. of codons."""
    which_orf = []
    for orf_id in orf:
        orf_bounds = orf[orf_id]
        orf_positions = sum([range(start,stop+1) for start,stop in orf_bounds], [])
        # pos residence in ORF or not in ORF?
        if loci in orf_positions:
            # relative loci, simply an index in the list of orf_positions ...
            relative_pos = orf_positions.index(loci)
            # the index of the codon it belongs to ...
            codon_idx = relative_pos//3
            # relative codon coordinates (indexes within CDS):
            relative_codon_coord = [codon_idx*3+i for i in range(3)]
            # loci is 1st,2nd or 3rd position in the codon?
            codon_shift = relative_pos - relative_codon_coord[0]
            # absolute codon coordinates (indexes in the genome):
            codon_coord = [orf_positions[_] for _ in relative_codon_coord]
            which_orf.append( (orf_id, (codon_idx,codon_shift,codon_coord) ) )
    return dict(which_orf)



# let's gather overall statistics ...
codon_position_mutated = []
codon_position_mutated_weight = []

dKD = []
dKD_weight = []



# as all of the alignments are loaded ...
# come up with some criteria for loci of interest ...
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
    # number of variable positions ...
    pos_with_var = descr.loc['freq']<descr.loc['count']
    print "There are %d variable positions to investigate for segment %s"%(pos_with_var.sum(),seg)
    # position indeces ...
    pos_ids = aln.loc[:,pos_with_var].columns
    # ORF of this segment ...
    seg_orf = ph1n1.orf[seg+'_pH1N1']
    seg_seq = ph1n1.genome[seg+'_pH1N1']
    #
    for pos in pos_ids:
        # unpack codon information ...
        prod_codon = get_loci_orf_assoc(seg_orf,pos)
        for product in prod_codon:
            codon_idx, codon_shift, codon_coords = prod_codon[product]
            # the consensus codon, AA and KD ...
            codon_itself = ''.join(seg_seq[i] for i in codon_coords)
            aa_itself = genetic_code[codon_itself]
            KD_itself = KD[aa_itself]
            # SNPs at this position ...
            nts_present = aln.loc[:,pos].value_counts().to_dict()
            ############################################
            # CODON SHIFT STATISTICS GATHERING ...
            codon_position_mutated.append(codon_shift)
            # non-consensus percent, same as 1.0-descr.loc['freq_ratio',pos]
            weight = sum(sorted(nts_present.values())[:-1])*1.0/sum(nts_present.values())
            assert abs(weight - (1.0-descr.loc['freq_ratio',pos]))<0.00000001
            codon_position_mutated_weight.append(weight)
            # CODON SHIFT STATISTICS GATHERING ...
            ############################################
            #
            # hence, possible codons are:
            possible_codons = []
            for snp in nts_present:
                mut_codon = list(codon_itself)
                mut_codon[codon_shift] = snp
                possible_codons.append(mut_codon)
                # STATISTICS GATHERING ...
                the_alt_codon = ''.join(mut_codon)
                if the_alt_codon != codon_itself:
                    the_alt_aa = genetic_code[the_alt_codon]
                    the_alt_KD = KD[the_alt_aa]
                    weight     = nts_present[snp]*1.0/sum(nts_present.values())
                    #
                    dKD.append(the_alt_KD - KD_itself)
                    dKD_weight.append(weight)
            # # amino acids ...
            # print "Product %s, position %d in protein (codon %s, aa %s)"%(product, codon_shift+1, codon_itself, genetic_code[codon_itself])
            # other_possible_codons = [''.join(codon) for codon in possible_codons if ''.join(codon)!=codon_itself]
            # print "outcome AAs are: %s"%str([genetic_code[ccc] for ccc in other_possible_codons]).strip('[]').replace('\'','')
            # # print str([genetic_code[ccc] for ccc in other_possible_codons]).strip('[]').replace('\'','')
            # print "dKD for AA subs: %s"%str([ KD[genetic_code[ccc]]-KD[genetic_code[codon_itself]] for ccc in other_possible_codons]).strip('[]')



        







































