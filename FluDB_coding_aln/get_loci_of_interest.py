# import re
import os
import sys
import copy
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
out_fname = 'test_loci.txt'
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





seg_aln_info = {}
# as all of the alignments are loaded ...
# come up with some criteria for loci of interest ...
# let's get some brief statistical info at first ...
for seg in sorted(segment_aln):
    aln = np.array([list(rec) for rec in segment_aln[seg]], np.character)
    aln = pd.DataFrame(aln)
    # replace gaps with None, so pd.isnull is the way to check gaps ...
    aln = aln.where(aln!='-',None)
    aln_descr = aln.describe().transpose()
    aln_descr['freq_ratio'] = aln_descr['freq']/aln_descr['count'] # freq_ratio ignores gaps ...
    aln_descr['gaps'] = aln.isnull().sum()
    aln_descr['variation'] = aln_descr['freq']<aln_descr['count']
    print "%s: %d variable positions"%(seg, aln_descr['variation'].sum())
    # position indeces ...
    loci_variable = aln_descr[aln_descr['variation']].index
    # ORF of this segment ...
    ref_seg_orf = ph1n1.orf[seg+'_pH1N1']
    ref_seg_seq = ph1n1.genome[seg+'_pH1N1']
    # detailed SNP info storage ...
    detailed_snp_info = {}
    detailed_snp_info['genome_pos'] = []
    detailed_snp_info['nt_ref'] = []
    detailed_snp_info['nt_snp'] = []
    detailed_snp_info['snp_freq'] = []
    detailed_snp_info['product'] = []
    detailed_snp_info['codon_shift'] = []
    detailed_snp_info['codon_ref'] = []
    detailed_snp_info['codon_snp'] = []
    detailed_snp_info['codon_pos'] = []
    # go over all genomic loci to gather SNPs information ...
    for pos in loci_variable:
        # unpack codon information ...
        prod_codon = get_loci_orf_assoc(ref_seg_orf,pos)
        for product in prod_codon:
            codon_idx, codon_shift, codon_coords = prod_codon[product]
            # the reference codon (as a list, first)...
            codon_ref = [ref_seg_seq[i] for i in codon_coords]
            # reference nucleotide and SNPs at this position ...
            nt_ref = codon_ref[codon_shift] # = aln_descr['top'][pos] = ref_seg_seq[pos]
            nts_present = aln.loc[:,pos].value_counts().to_dict()
            # hence, possible codons are:
            for nt_snp in nts_present:
                if nt_snp != nt_ref:
                    codon_snp = copy.deepcopy(codon_ref)
                    codon_snp[codon_shift] = nt_snp
                    # at the innermost level let's gather detailed SNP info here ...
                    detailed_snp_info['genome_pos'].append(pos)
                    detailed_snp_info['nt_ref'].append(nt_ref)
                    detailed_snp_info['nt_snp'].append(nt_snp)
                    detailed_snp_info['snp_freq'].append(nts_present[nt_snp])
                    detailed_snp_info['product'].append(product)
                    detailed_snp_info['codon_shift'].append(codon_shift)
                    detailed_snp_info['codon_ref'].append(''.join(codon_ref))
                    detailed_snp_info['codon_snp'].append(''.join(codon_snp))
                    detailed_snp_info['codon_pos'].append(codon_idx)
    ###################################################################################
    detailed_snp_info = pd.DataFrame(detailed_snp_info)
    # SNPs information is gathered ...
    # let's merge 'aln_descr' and 'detailed_snp_info' dataframes
    seg_aln_info[seg] = pd.merge(aln_descr,detailed_snp_info,left_index=True,right_on='genome_pos',how='outer').reset_index(drop=True)




for seg in seg_aln_info:
    # aln_info is merely a reference to object
    # so, changing stuff by ref, will change the object itself ...
    seg_aln_info[seg] = seg_aln_info[seg].reset_index(drop=True)
    aln_info = seg_aln_info[seg]
    aln_info['genome_pos'] = aln_info['genome_pos'].astype('int')
    aln_info['seg'] = seg
    # reconstruct AAs by the codons ...
    aln_info['aa_ref'] = aln_info['codon_ref'].map(genetic_code)
    aln_info['aa_snp'] = aln_info['codon_snp'].map(genetic_code)
    # get KDs by AAs ...
    aln_info['KD_ref'] = aln_info['aa_ref'].map(KD)
    aln_info['KD_snp'] = aln_info['aa_snp'].map(KD)
    # delta KD, snp minus reference ...
    aln_info['dKD'] = aln_info['KD_snp'] - aln_info['KD_ref']


result_df = pd.concat( [seg_aln_info[seg][seg_aln_info[seg]['variation']][['seg','genome_pos']] for seg in sorted(seg_aln_info)], ignore_index=True )
result_df = result_df.drop_duplicates().reset_index(drop=True)
# then we would need some criteria to select interesting loci and just output them ...
# let's say now we're interested in all of the variable loci:
result_df.to_csv(out_fname,index=False)




































