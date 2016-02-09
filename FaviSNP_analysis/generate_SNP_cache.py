import sys
import re
import os
# import string
import numpy as np
# import pandas as pd
# import subprocess as sub
from Bio import SeqIO
from Bio.Data import CodonTable

# import scipy.sparse as sparse
# import scipy.stats as st

# import flu_module as flu




class influenza(object):
    """influenza is a class describing all possible feature of a given influenza genome,
    including genome itself, lengths of segments, ORFs, RNA profiles etc."""
    #
    #
    #
    def __init__(self, genome_fname, ORF_fname):
        # read influenza genome ...
        self.genome = dict([(seq.id,seq.seq) for seq in SeqIO.parse(genome_fname,"fasta")])
        # length of the segments ...
        self.seglen = dict([(segid,len(self.genome[segid])) for segid in self.genome])
        # genome length ...
        self.genlen = sum(self.seglen.values())
        # number of segments ...
        self.segnum = len(self.genome)
        # name of the file with the genome ...
        self.genome_fname = genome_fname
        # just segment's names/ids - whatever ...
        self.segments = sorted(self.seglen.keys())
        # proteins ...
        # initialize ORF data ...
        ORF_per_segment = {}
        for seg in self.segments:
            ORF_per_segment[seg] = {}
        # read ORF data ...
        with open(ORF_fname,'r') as fp:
            for line in fp.readlines():
                line = line.strip().split()
                # bound must be specified with int numbers ...
                # shift them left by 1, to make them 0-numbered ...
                orf_bounds = [int(_)-1 for _ in line[2:]]
                # >0, even # of numbers ...
                if (len(orf_bounds)%2>0)or(len(orf_bounds)==0):
                    raise ValueError("ORF boundaries must be specified with even number of positions!")
                # form pairs: [1,2,3,4] -> [(1,2),(3,4)]
                orf_bounds = zip(orf_bounds[0::2],orf_bounds[1::2])
                seg = line[0]+'FULL' # for some reason, segments in .orf files are mentioned without FULL addition ...
                # the keys must coincide with the existing ones ...
                assert (seg in ORF_per_segment.keys())
                # store ORF's translation(protein) name with the boundaries ... 
                protein = line[1]
                ORF_per_segment[seg][protein] = orf_bounds
        self.orf = ORF_per_segment



def isin_orf(orf,pos):
    """checks if a given position is in any of the ORFs.
    Returns list of the ORFs that it participates in with the codon coords."""
    which_orf = []
    for orf_id in orf:
        orf_bounds = orf[orf_id]
        # assuming total length of the ORF is ~3
        assert sum((stop-start+1) for start,stop in orf_bounds)%3 == 0
        # if that holds, start analysis ...
        orf_positions = sum([range(start,stop+1) for start,stop in orf_bounds], [])
        # pos residence in ORF or not in ORF?
        if pos in orf_positions:
            # relative pos, simply an index in the list of orf_positions ...
            relative_pos = orf_positions.index(pos)
            # the index of the codon it belongs to ...
            codon_idx = relative_pos//3
            # relative codon coordinates (indexes within CDS):
            relative_codon_coord = [codon_idx*3+i for i in range(3)]
            # pos is 1st,2nd or 3rd position in the codon?
            codon_shift = relative_pos - relative_codon_coord[0]
            # absolute codon coordinates (indexes in the genome):
            codon_coord = [orf_positions[_] for _ in relative_codon_coord]
            which_orf.append( (orf_id, (codon_idx,codon_shift,codon_coord) ) )
    return dict(which_orf)



# def isin_orf(orf,pos):
#     """checks if a given position is in any of the ORFs.
#     Returns list of the ORFs that it participates in with the codon coords."""
#     which_orf = []
#     for orf_id in orf:
#         orf_bounds = orf[orf_id]
#         # assuming total length of the ORF is ~3
#         assert sum((stop-start+1) for start,stop in orf_bounds)%3 == 0
#         # if that holds, start analysis ...
#         for start,stop in orf_bounds:
#             # NOT CORRECT TO BE EDITED ................!!!!!!!!!!!!!!!!!!!
#             if start <= pos <= stop:
#                 # which codon ...
#                 codon = (pos-start)//3
#                 # codon coordinates are ...
#                 codon_start = codon*3 + start
#                 # c1 = codon*3 + start + 1
#                 codon_stop  = codon*3 + start + 2
#                 # which nucleotide position "pos" corresponds to?
#                 codon_shift = pos - codon_start
#                 # store stuff in  ...
#                 which_orf.append( (orf_id, (codon,codon_shift,slice(codon_start,codon_stop+1)) ) )
#     return dict(which_orf)



# genome_fname = "../FluGenomes/H1N1FULL.fa"
genome_fname = "../FluGenomes/H1N1FULL.fa"
orf_fname = "../FluGenomes/xlator_data/H1N1FULL.orf"
#


h1n1 = influenza(genome_fname,orf_fname)

genetic_code = CodonTable.standard_dna_table.forward_table
stop_codons = dict([ (codon,'*') for codon in CodonTable.standard_dna_table.stop_codons ])
# make genetic code contain both 61 aacid codons and 3 stop codons ...
genetic_code.update(stop_codons)

# mutate genome position
for seg in h1n1.segments:
    seg_len = h1n1.seglen[seg]
    segment = h1n1.genome[seg]
    the_orf = h1n1.orf[seg]
    for pos in xrange(seg_len):
        ref_nt  = segment[pos]
        mut_nucs = set('ATGC')-set(ref_nt)
        inorf = isin_orf(the_orf,pos)
        if not inorf:
            for mut_nt in mut_nucs:
                # example: seg1_H1N1FULL_A456G NORF
                line_out = "%s_%s%d%s NORF"%(seg,ref_nt,pos+1,mut_nt)
                print line_out
        else:
            # might be inefficient, but works for our output ...
            for mut_nt in mut_nucs:
                line_out = "%s_%s%d%s"%(seg,ref_nt,pos+1,mut_nt)
                for product in inorf:
                    aa_idx, codon_shift, codon_coord = inorf[product]
                    # we use explicit coordinates for all 3 nucs in codon,
                    # to handle within codon splice sites properly ...
                    the_codon = ''.join(segment[_] for _ in codon_coord)
                    ref_aa = genetic_code[the_codon]
                    # modifying the ref codon to mutant ...
                    mut_codon = list(the_codon)
                    mut_codon[codon_shift] = mut_nt
                    mut_codon = ''.join(mut_codon)
                    mut_aa = genetic_code[mut_codon]
                    # print ...
                    line_out += " %s_%s%d%s"%(product,ref_aa,aa_idx+1,mut_aa)
                print line_out






























