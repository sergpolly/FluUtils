import sys
# import re
# import os
# import string
import random as rnd
import numpy as np
import pandas as pd
import subprocess as sub
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio import SeqRecord
# import Bio
import copy
import re

import scipy.sparse as sparse


import matplotlib as mpl
mpl.rcParams['lines.linewidth'] = 4
# from matplotlib.patches import Polygon
import matplotlib.transforms as transforms

#Font sizes and stuff ...
mpl.rcParams['font.family']=['FreeSans']
# fig,ax = plt.subplots(nrows=1,ncols=1,figsize=(9,6))
# # ax = fig.add_subplot(111)
# plt.subplots_adjust(left=0.1,right=0.80,bottom=0.1,top=0.9) 
fontsize_xy_label = 21
fontsize_xy_ticks_labels = 20
fontsize_useles_title = 15
linewidth_plot = 6


# before fixing the UTR problem seg lengthes were :
# seg_lengths = {'seg1_H1N1':2314, 'seg2_H1N1':2302, 'seg3_H1N1':2203, 'seg4_H1N1':1776, 'seg5_H1N1':1497, 'seg6_H1N1':1427, 'seg7_H1N1':1006, 'seg8_H1N1':870}



#####################################
#####################################
#####################################
class influenza(object):
    """influenza is a class describing all possible feature of a given influenza genome,
    including genome itself, lengths of segments, ORFs, RNA profiles etc."""
    #
    #
    #
    def __init__(self, genome_fname, ORF_fname, mutations=None):
        # read influenza genome ...
        self.genome = dict([ (seq.id,seq.seq) for seq in SeqIO.parse(genome_fname,"fasta") ])
        # length of the segments ...
        self.seglen = dict( [(segid,len(segseq)) for segid,segseq in self.genome.iteritems() ])
        # genome length ...
        self.genlen = sum(self.seglen.values())
        # number of segments ...
        self.segnum = len(self.genome)
        # name of the file with the genome ...
        self.genome_fname = genome_fname
        # just segment's names/ids - whatever ...
        self.segments = sorted(self.seglen.keys())
        # (-) sense sequence is needed for RNA folding ...
        self.genome_comp = dict([ (segid,segseq.complement()) for segid,segseq in self.genome.iteritems() ])
        # proteins ...
        # self.segcode = dict([(i,name) for i,name in enumerate(['PB2','PB1','PA','HA','NP','NA','M1/M2','NS1/NEP'])])
        # initialize ORF data using segment id-s from genome_fname file ...
        ORF_per_segment = {}
        for seg in self.segments:
            ORF_per_segment[seg] = {}
        # read ORF data ...
        with open(ORF_fname,'r') as fp:
            for line in fp.readlines():
                line = line.strip().split()
                #####################################################
                #  PERFORM SOME SIMPE ORF CHECKS ALONG ...
                #####################################################
                seg = line[0] # make sure you have the same segment ids in genome_fname and ORF_fname ...
                # the keys must coincide with the ones from genome_fname file ...
                if seg not in self.segments:
                    raise ValueError("segment ids in %s must coincide with the ids in %s!"%(ORF_fname,genome_fname))
                # bound must be specified with int numbers ...
                # shift them left by 1, to make them 0-numbered ...
                orf_bounds = [int(_)-1 for _ in line[2:]]
                # >0, even # of numbers ...
                if (len(orf_bounds)%2>0)or(len(orf_bounds)==0):
                    raise ValueError("ORF boundaries must be specified with even number of coordinates!")
                # form pairs: [1,2,3,4] -> [(1,2),(3,4)]
                orf_bounds = zip(orf_bounds[0::2],orf_bounds[1::2])
                # check 3-multiplicity for the ORF length ...
                if sum((oe-ob+1) for ob,oe in orf_bounds)%3 != 0:
                    raise ValueError("ORF length must be multiple of 3! Corrupted ORF prevents codons indentification!")
                # store ORF's translation(protein) name with the boundaries ... 
                protein_id = line[1]
                ORF_per_segment[seg][protein_id] = orf_bounds
        self.orf = ORF_per_segment
    #
    #
    #
    # def __RNAfold(self,seq_dict,maxBPspan=-1):
        # inner function to turn dict of sequences to fasta file (or string) ...
        def local_seq_dict_fasta(seq_dict):
            return '\n'.join( ">%s\n%s"%(sid,seq) for sid,seq in seq_dict.iteritems() )
        # inner function to parse output of RNAfold ...
        def local_parse_RNAfold(output_str):
            mfe_pat = re.compile("\([^()]+\d+\.\d+\)") # stuff like '( -9.10)' or '(  0.10)' ...
            out_list = output_str.strip().split('\n') # list of ~3 size: id, seq, RNA structure with MFE
            sid_struct_mfe_list = [] # storage for results ...
            # then loop in three-s:
            # zip(l[::3],l[1::3],l[2::3]) turns l=[1,2,3,4,5,6] to [(1,2,3),(4,5,6)] ...
            for sid,seq,struct in zip( out_list[::3],out_list[1::3],out_list[2::3] ):
                sid = sid.strip().strip('<>')
                search_res = mfe_pat.search(struct)
                mfe = float(search_res.group().strip('()'))
                struct = struct[:search_res.start()].strip()
                sid_struct_mfe_list.append( (sid,(struct,mfe)) )
            return dict(sid_struct_mfe_list)
        # now launch RNAfold and collect results ...
        cmd = ['RNAfold','--maxBPspan',str(maxBPspan)]
        RNAfold_process = sub.Popen(cmd,stdin=sub.PIPE,stdout=sub.PIPE)
        output,errput = RNAfold_process.communicate(local_seq_dict_fasta(seq_dict))
        struct_mfe_dict = local_parse_RNAfold(output)
        # returning parsed structures with MFE ...
        return struct_mfe_dict
    #
    #
    #
    def __RNALfold(self,seq_dict,window=30):
        # inner function to turn dict of sequences to fasta file (or string) ...
        def local_seq_dict_fasta(seq_dict):
            return '\n'.join( ">%s\n%s"%(sid,seq) for sid,seq in seq_dict.iteritems() )
        # inner function to parse output of RNAfold (tested & working) ...
        def local_parse_RNALfold(output_str):
            hairpin_pat = re.compile("^([\.\(\)]+)\s+\(([^()]+\d+\.\d+)\)\s+(\d+)$") # stuff like ".(((((......))))). ( -6.00)  25" ...
            out_list_struct_seq = output_str.strip('<>').split('>') # divide them by sequences first ...
            sid_list = [] # outer storage for results ...
            for local_struct_str in out_list_struct_seq:
                hairpin_list = local_struct_str.strip().split('\n')
                sid = hairpin_list[0].strip()
                mfe_total = float(hairpin_list[-1].strip().strip('()'))
                seq = hairpin_list[-2]
                # the rest: hairpin_list[1:-2] are local structures (hairpins with their mfe and positions) ...
                parsed_hairpins = []
                for hairpin_descr in hairpin_list[1:-2]:
                    hairpin_struct,hairpin_mfe,hairpin_pos = hairpin_pat.match(hairpin_descr).groups()
                    hairpin_mfe = float(hairpin_mfe)
                    hairpin_pos = int(hairpin_pos)
                    parsed_hairpins.append( (hairpin_pos,hairpin_struct,hairpin_mfe,mfe_total) )
                sid_list.append( (sid,parsed_hairpins) )
            return dict(sid_list)
        # now launch RNALfold and collect results ...
        cmd = ['RNALfold','-L',str(window)]
        RNALfold_process = sub.Popen(cmd,stdin=sub.PIPE,stdout=sub.PIPE)
        output,errput = RNALfold_process.communicate(local_seq_dict_fasta(seq_dict))
        local_struct_dict = local_parse_RNALfold(output)
        # returning parsed structures with MFE ...
        return local_struct_dict
    #
    #
    # Local MFE example output ...
    ###################################################################
    # >seg8_H1N1FULL
    # .(((((((((...)))))..)))). ( -7.90)  123
    # .((((.....(((((...))))).....)))). (-10.30)  118
    # .((....)). ( -1.60)    1
    # AGCGAAAGCAGGGUGGCAAAGACAUAAUGGAUUCCCACACUGUGUAC...
    #  (-199.00)
    ###################################################################
    # segments are overlapping or not ?!
    def __get_overlap(self,interval1,interval2):
        idx1,start1,str_len1 = interval1
        idx2,start2,str_len2 = interval2
        # L = max(A.left,B.left)
        # R = min(A.right,B.right)
        # L <= R - <=> intersection!!! overlap
        L = max(start1, start2)
        R = min(start1+(str_len1-1), start2+(str_len2-1))
        return (L <= R)
    #
    def __get_overlap_matrix(self,hairpin_df):
        # this function tested and produces meaningful results ...
        # BEWARE!
        # we don't need overlap matrix itself other than for illustration purposes
        def local_get_overlap(starts,lens):
            # check universality of the function ...
            start1,start2 = starts
            len1,len2 = lens
            # L = max(A.left,B.left)
            # R = min(A.right,B.right)
            # L <= R - <=> intersection!!! overlap
            L = max(start1, start2)
            R = min(start1+(len1-1), start2+(len2-1))
            return (L <= R)
        #
        hairpin_list = list( hairpin_df[['start','struct_len']].itertuples() )
        hairpin_list_len = len(hairpin_list)
        #
        # 'hairpin_list_len' by 'hairpin_list_len' matrix of zeros to be filled in ...
        overmat = np.zeros((hairpin_list_len,hairpin_list_len),dtype=np.int)
        for i in range(hairpin_list_len):
            for j in range(i+1,hairpin_list_len):
                idx,istart,ilen = hairpin_list[i]
                jdx,jstart,jlen = hairpin_list[j]
                # if these hairpins are overlaping ...
                if local_get_overlap( (istart,jstart), (ilen,jlen) ):
                    overmat[idx,idx] = 1
                    overmat[jdx,jdx] = 1
                    overmat[idx,jdx] = 1
                    overmat[jdx,idx] = 1
        return overmat
    #
    #
    def __get_constraint_matrix(self,hairpin_df):
        # this function tested and produces meaningful results ...
        # BEWARE!
        # Constraint matrix is different from symmetric overlap matrix.
        # Constraint one is needed to solve MIP (optimization) problem.
        # see github wiki for elaborate discussion ...
        def local_get_overlap(starts,lens):
            # check universality of the function ...
            start1,start2 = starts
            len1,len2 = lens
            # L = max(A.left,B.left)
            # R = min(A.right,B.right)
            # L <= R - <=> intersection!!! overlap
            L = max(start1, start2)
            R = min(start1+(len1-1), start2+(len2-1))
            return (L <= R)
        #
        hairpin_list = list( hairpin_df[['start','struct_len']].itertuples() )
        hairpin_list_len = len(hairpin_list)
        #
        constramat = []
        for i in range(hairpin_list_len):
            for j in range(i+1,hairpin_list_len):
                idx,istart,ilen = hairpin_list[i]
                jdx,jstart,jlen = hairpin_list[j]
                # if these hairpins are overlaping ...
                if local_get_overlap( (istart,jstart), (ilen,jlen) ):
                    # THIS IS VERY IMPORTANT: PROPER WAY TO FILL IN CONSTRAINT MATRIX  ...
                    # embrace pair-wise constraint description. Simply list all mutatually exclusive hairpins ...
                    tmp_vec = np.zeros(hairpin_list_len,dtype=np.int)
                    tmp_vec[idx] = 1
                    tmp_vec[jdx] = 1
                    constramat.append(tmp_vec)
        # to matrix from a list of vectors ...
        constramat = np.asarray(constramat,dtype=np.int)
        # dense matrix is filled. Get the COO representation for the matrix ...
        constramat_coo = sparse.coo_matrix(constramat)
        # 'row': constramat_coo.row
        # 'col': constramat_coo.col
        # 'val': constramat_coo.data
        # pairwise overlaps are established, returing COO matrix ...
        return constramat_coo
    #
    #
    #
    def __optimize_overlap_glpk(self,overmat,objective,uniq_id,bip_path='./rnaglpk/svbin'):
        fname_obj = "%s_objective.dat"%uniq_id
        fname_cst = "%s_sparse.dat"%uniq_id
        # number of constraints ...
        overmat_size = len(overmat.data)
        objective_size = objective.size
        ##################################################
        # writing the sparse matrix of the constraints ...
        with open(fname_cst,'w') as fp:
            # let's write the total number of lines at first:
            fp.write("%d\n"%overmat_size)
            # sparse matrix of constraints is to follow (mind 1-based indexing of GLPK) ...
            for index,i,j,v in zip(range(overmat_size),overmat.row, overmat.col, overmat.data):
                fp.write("%d a_%d_%d %d %d %d\n" % (index+1,i+1,j+1,i+1,j+1,v))
        ###################################################
        # objective coefficients ...
        with open(fname_obj,'w') as fp:
            # total number of items first ...
            fp.write("%d\n"%objective_size)
            # all items second ...
            for index,obj_coeff in enumerate(objective):
                fp.write("%d a_%d %.2f\n" % (index+1,index+1,obj_coeff))
        # run the optimizer ...
        cmd = "%s %s %s"%( os.path.join(bip_path,'run_bip'), fname_obj, fname_cst )
        result = sub.check_output(cmd,shell=True)
        # read last 'objective_size' lines from the result ...
        result_out = np.asarray( [ int(i) for i in result.strip().split('\n')[-objective_size:] ] )
        # return included hairpir coefficients ...
        return result_out
    #
    #
    #
    #
    def __RNALfold_opt(self,seq_dict,dg_NP,windows=[24,36,]):
        # calculate local RNA structure using RNALfold ...
        df_mfe = []
        for length in windows:
            # tmp_mfe = self.__RNALfold(seq_dict,window=length)
            tmp_mfe = __RNALfold(self,seq_dict,window=length)
            for sid in seq_dict:
                tmp_df = pd.DataFrame(tmp_mfe[sid],columns=['start','struct','mfe','opt_mfe'])
                tmp_df['window'] = length
                tmp_df['seq_id'] = sid
                df_mfe.append(tmp_df)
        lfold_mfe = pd.concat(df_mfe,ignore_index=True)
        ###################################################
        # WARNING! some structures include unpaired flanks .(). , some of them do not: .() or (). 
        # all that affects start position and struct_len.
        #
        # (1) adjust starting position accordingly: +1 for all .() and .().; keep 1-indexing for now...
        lfold_mfe['start'] += lfold_mfe['struct'].str.startswith('.').apply(int)
        # (2) strip those unpaired flanks '.', to enforce matching ends for all ((...))
        lfold_mfe['struct'] = lfold_mfe['struct'].str.strip('.')
        # (3) get the structure length finally ...
        lfold_mfe['struct_len'] = lfold_mfe['struct'].str.len()
        # (4) calculate hairpin MFE or dG per nucleotide ...
        lfold_mfe['mfe_per_nt'] = lfold_mfe['mfe']/lfold_mfe['struct_len']
        # (5) calculate objective coefficients for optimizaiton (MFE corrected to NP attraction)...
        lfold_mfe['objective_mfe'] = (lfold_mfe['mfe_per_nt'] - dg_NP)*lfold_mfe['struct_len']
        # lfold_mfe is a mixed DataFrame with all the seq_id and all window lengths ...
        #
        #
        # ... SEVERAL SCENARIOS ARE POSSIBLE FROM THIS POINT ...
        # (1) drop all hairpoins yielded by Lfolding with different window lengths (Do this before, probably) ...
        lfold_mfe = lfold_mfe.drop_duplicates(subset=[col for col in lfold_mfe.columns if (col not in ['window','opt_mfe']) ])
        # (2) drop those hairpins, melted by NP ...
        lfold_mfe = lfold_mfe[ lfold_mfe['objective_mfe']<0.0 ]
        # (3) reset index, just in case ...
        lfold_mfe = lfold_mfe.reset_index(drop=True)
        #
        # Next, we would need to fill in overlap matrix ...
        lfold_sid = lfold_mfe.groupby('seq_id')
        #
        # UNDER CONSTRUCTION ....
        return lfold_sid
        #
        #
        # loop over structure/set of hairpins corresponding to a particular 'sid'
        for sid in lfold_sid.groups:
            lfold_mfe_sid = lfold_sid.get_group(sid).sort(columns='start').reset_index(drop=True)
            overmat_sid = self.__get_constraint_matrix(lfold_mfe_sid) # modify later, to avoid copying the whole 'lfold_mfe_sid' ...
            # overmat is a dict with 'row','col' and 'val' keys ...
            hairpin_opt_list = self.__optimize_overlap_glpk(overmat_sid,lfold_mfe_sid['objective_mfe'],sid,bip_path='./rnaglpk/svbin') # change path later on ...
            # 0,1-mask vector describing hairpins in optimal structure ...
        #
        ################################################################################
        #   OLD STUFF FOLLOWS, TO BE REWRITTEN ...
        # path = r"./Dropbox (UMASS MED - BIB)/CellSupProject/cellsup2015_with_RNA/rnaglpk/svbin"
        # path = r"./Dropbox\ \(UMASS\ MED\ -\ BIB\)/CellSupProject/cellsup2015_with_RNA/rnaglpk/svbin"
        ################################################################################
        # now store that MFE object dictionary as an influenza class attribute ...
        self.optimal_MFE = mfe
        #
        #
        # we'd also calculate the dg profile over each segment: 0 if NP-bound, dg(per nucleotide) if belongs to some stable structure ...
        #
        dg_profile = {}
        pairing_profile = {}
        loops_profile = {}
        mfe_loops = {}
        #
        pairing_dict = {'.':1,'(':2,')':2}
        # big hairpin loops of 4 nt and more ...
        big_loop = re.compile('\.{4,}')
        #
        # now let's try to plot that stuff ...
        for seg in self.segments:
            seglen = self.seglen[seg]
            # take only those RNA that've been selected by optimizer ...
            seg_struct = mfe[seg][mfe[seg]['inc_rna']==1]
            # let's generate RNA profiles based on the optimized results ...
            # print seglen
            dg_profile[seg] = np.zeros(seglen)
            pairing_profile[seg] = np.zeros(seglen,dtype = np.int)
            loops_profile[seg] = np.zeros(seglen,dtype = np.int)
            mfe_loops[seg] = []
            #
            for dat in seg_struct.iterrows():
                struct = dat[1]['struct']
                start = dat[1]['start']
                str_len = dat[1]['str_len']
                dg_nuc = dat[1]['dg_nuc']
                # search big loops ...
                for res in big_loop.finditer(struct):
                    # the profile of big loops, those potentially available for intersegmental pairing ...
                    loops_profile[seg][start+res.start():start+res.end()] = 1
                    #  sequences of those big loops, just to look at who they really are ...
                    mfe_loops[seg].append(( start+res.start(), self.genome[seg][start+res.start():start+res.end()] ))
                    #
                dg_profile[seg][start:start+str_len] = dg_nuc
                pairing_profile[seg][start:start+str_len] = [ pairing_dict[char] for char in struct[1:-1] ]
                # maybe struct[1:-1] if 0 and -1 are '.'
        #
        # now store that as an influenza class attribute ...
        self.dg_profile = dg_profile
        self.pair_profile = pairing_profile
        self.loops_profile = loops_profile
        self.loops = mfe_loops





    #
    #
    #
    #
    # calculate local RNA profile of a flu genome taking NP attraction into account(dg_NP threshold) ... 
    def calculate_rna_lMFE_optimized(self, dg_NP, windows=[24,36,]):
        #
        genome_fname=self.genome_fname
        #
        mfe = {}
        for seg in self.segments:
            mfe[seg] = {}
        #
        # calculate local MFE structures for several window sizes,
        # to get a good grip on possible stable structures and their substructures ...
        for length in windows:
            self.calculate_rna_lMFE(length=length)
            for seg in self.segments:
                mfe[seg][length] = self.local_MFE[seg]
        # structures are stored in 'mfe':
        # organized first by segment, next by the window size
        #
        # next, we are calculating an overlap matrix, and store it in COO format:
        overlap_mat = {}
        overlat_sparse_mat = {}
        for seg in self.segments:
            #
            # mixing different folding lengths to find optimal folding solution ...
            mfe[seg] = pd.DataFrame(
                [(v1,v2,v3,v4,key) for key,values in mfe[seg].iteritems() for v1,v2,v3,v4 in values],
                columns = ['struct','start','str_len','dG','lfold'])
            #
            # sort by the starting point ... (1-indexed)
            mfe[seg] = mfe[seg].sort(columns=['start',]).reset_index(drop=True)
            # IMPORTANT adjustments to the segment dataframe ...
            # each structured rna segment includes 2 unpaired nucleotides (left & right)
            # we adjust corresponding segment lengths ...
            mfe[seg].str_len -= 2
            #
            # let's add energy per nt here
            mfe[seg]['dg_nuc'] = np.true_divide(mfe[seg]['dG'],mfe[seg]['str_len'])
            #
            # and objective coefficients along the way:
            # this is simply the dG of the structural element, corrected for the NP attraction ...
            mfe[seg]['objective_coeff'] = ( mfe[seg]['dg_nuc'] - dg_NP ) * mfe[seg]['str_len']
            #
            ##########################################################
            # clean redundant data:
            # drop duplicates first, ('lfold' is not taken inot account here ...)
            mfe[seg] = mfe[seg].drop_duplicates(subset=['struct', 'start', 'str_len', 'dG', 'dg_nuc', 'objective_coeff'])
            # drop those that are melted by NP,
            mfe[seg] = mfe[seg][mfe[seg]['objective_coeff']<0.0]
            # reset index - just in case ...
            mfe[seg] = mfe[seg].reset_index(drop=True)
            # looks like it's the best we can do so far...
            ##########################################################
            #
            # start indexes are actually OK, if you reconsider RNA sequences as 0-based indexed.
            #
            rna_seg_num = mfe[seg].shape[0]
            rna_str_seg = list(mfe[seg][['start','str_len']].itertuples())
            # integer matrix of overlap information.
            overlap_mat[seg] = []
            # overlap_mat[seg] = np.zeros((rna_seg_num,rna_seg_num),dtype=np.int)
            #
            for i in range(rna_seg_num):
                for j in range(i+1,rna_seg_num):
                    if self.__get_overlap(rna_str_seg[i],rna_str_seg[j]):
                        tmp_vec = np.zeros(rna_seg_num,dtype=np.int)
                        tmp_vec[i] = 1
                        tmp_vec[j] = 1
                        # overlap_mat[seg][i,i] = 1 # little bit redundant, but OK ...
                        # overlap_mat[seg][i,j] = 1
                        overlap_mat[seg].append(tmp_vec)
            overlap_mat[seg] = np.asarray(overlap_mat[seg],dtype=np.int)
            #
            overlat_sparse_mat[seg] = sparse.coo_matrix(overlap_mat[seg])
        ################################################################################
        ################################################################################
        ################################################################################
        # cycle through all of the segments ...
        for seg in self.segments:
            print "working on segment %s"%seg
            #
            fname_obj = "objective.dat"
            fname_cst = "sparse.dat"
            #
            # writing the sparse matrix of the constraints ...
            with open(fname_cst,'w') as fp:
                # glpk uses 1-based indexing ...
                # let's write the total number of lines at first:
                fp.write("%d\n" % len(overlat_sparse_mat[seg].data))
                #
                index = 1
                for i,j,v in zip(overlat_sparse_mat[seg].row+1,overlat_sparse_mat[seg].col+1,overlat_sparse_mat[seg].data):
                    fp.write("%d a%d%d %d %d %d\n" % (index,i,j,i,j,v))
                    index += 1
            ##############################
            with open(fname_obj,'w') as fp:
                # total number of items first ...
                fp.write("%d\n" % mfe[seg]['objective_coeff'].shape[0])
                # all items second ...
                for index,dg in mfe[seg]['objective_coeff'].iteritems():
                    fp.write("%d a%d %.2f\n" % (index+1,index+1,dg))
            # run the optimizer ...
            cmd_per_seg = "./rnaglpk/svbin/run_bip %s %s"%(fname_obj,fname_cst)
            result = sub.check_output(cmd_per_seg,shell=True)
            # ...
            expected_lines_useful = mfe[seg]['objective_coeff'].shape[0]
            result_out = np.asarray( [ int(i) for i in result.strip().split('\n')[-expected_lines_useful:] ] )
            # ...
            mfe[seg]['inc_rna'] = result_out
        #
        # now store that MFE object dictionary as an influenza class attribute ...
        self.optimal_MFE = mfe
        #
        #
        # we'd also calculate the dg profile over each segment: 0 if NP-bound, dg(per nucleotide) if belongs to some stable structure ...
        #
        dg_profile = {}
        pairing_profile = {}
        loops_profile = {}
        mfe_loops = {}
        #
        pairing_dict = {'.':1,'(':2,')':2}
        # big hairpin loops of 4 nt and more ...
        big_loop = re.compile('\.{4,}')
        #
        # now let's try to plot that stuff ...
        for seg in self.segments:
            seglen = self.seglen[seg]
            # take only those RNA that've been selected by optimizer ...
            seg_struct = mfe[seg][mfe[seg]['inc_rna']==1]
            # let's generate RNA profiles based on the optimized results ...
            # print seglen
            dg_profile[seg] = np.zeros(seglen)
            pairing_profile[seg] = np.zeros(seglen,dtype = np.int)
            loops_profile[seg] = np.zeros(seglen,dtype = np.int)
            mfe_loops[seg] = []
            #
            for dat in seg_struct.iterrows():
                struct = dat[1]['struct']
                start = dat[1]['start']
                str_len = dat[1]['str_len']
                dg_nuc = dat[1]['dg_nuc']
                # search big loops ...
                for res in big_loop.finditer(struct):
                    # the profile of big loops, those potentially available for intersegmental pairing ...
                    loops_profile[seg][start+res.start():start+res.end()] = 1
                    #  sequences of those big loops, just to look at who they really are ...
                    mfe_loops[seg].append(( start+res.start(), self.genome[seg][start+res.start():start+res.end()] ))
                    #
                dg_profile[seg][start:start+str_len] = dg_nuc
                pairing_profile[seg][start:start+str_len] = [ pairing_dict[char] for char in struct[1:-1] ]
                # maybe struct[1:-1] if 0 and -1 are '.'
        #
        # now store that as an influenza class attribute ...
        self.dg_profile = dg_profile
        self.pair_profile = pairing_profile
        self.loops_profile = loops_profile
        self.loops = mfe_loops
    #
    #
    #
    #
    #
    #
    #
    #
    #
    # 
    # BROWSER WILL BREAK MOST LIKELY, BECAUSE ORF's DETERMINATION HAS BEEN UPDATED ...
    def genome_browser(self):
        """A function returning figure and axis for the influenza genome browser.
        It generates 8 empty panels of lengths proportional to the segment and
        allows """
        fig = plt.figure(figsize=(12,12))
        #
        # axs = [None for _ in range(self.segnum)]
        #
        # horizontal ...
        width_longestseq = 0.9
        left = (1. - width_longestseq)/2.
        # vertical ...
        bottom_lowest = 0.05
        v_delta = bottom_lowest
        span_to_height_coeff = 0.8
        num_of_spans = self.segnum - 1
        height = (1.-2.*v_delta)/(self.segnum + span_to_height_coeff*num_of_spans)
        v_span = span_to_height_coeff*height
        #
        #
        bottom = bottom_lowest
        max_seglen = float(max(self.seglen.values()))
        axs = {}
        # starting from the longest segment down to the shortest one ...
        for seg,seglen in sorted(self.seglen.items(),key=lambda x: (-x[1],x[0]),reverse=True):
            coeff = seglen/max_seglen
            axs[seg] = fig.add_axes([left, bottom, width_longestseq*coeff, height])
            bottom += height + v_span
        #
        #
        #go through axes to setup all the limits and  stuff ...
        # for i,ax in enumerate(axs):
        for seg,ax in axs.items():
            # DRAW ORF AT EVERY AXIS ...
            a,b = self.orf[seg]
            ax.axvspan(xmin=a,xmax=b,alpha=0.9, facecolor='0.9', edgecolor='none',zorder=1)
            #
            # draw middle lines ...
            ax.plot([0,1],[0.5,0.5],linestyle='-',color='k',linewidth=1.,transform=ax.transAxes)
            # #
            # # fontsize = 8
            x_step = 500 if self.seglen[seg]>=1500 else 200
            ax.set_xlim((0, self.seglen[seg] ))
            # no axis labels ...
            # ax.set_xlabel('nucleotide',fontsize=fontsize_xy_label-3)
            # ax.set_ylabel('SNP group',fontsize=fontsize_xy_label-3)\
            trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
            ax.xaxis.set_label_coords(500,-0.35*0.5,transform=trans)
            ax.set_xticks(np.arange(0, self.seglen[seg] , x_step))
            ax.set_xticklabels([1]+range(x_step,self.seglen[seg],x_step),fontsize=fontsize_xy_ticks_labels)
            # #
            # #
            # # we dont know much about y-range and limits - it'll be controlled outside!
            # # it'll depend on the type of data that we'll plot there ...
            # ylim = 1
            # ax.set_ylim((-1.*ylim,ylim))
            # yticks_range = range(ylim+1)
            # # clean y-axis just for clarity ...
            ax.set_yticks([])
            ax.set_yticklabels([])
            # # the segment length bar ...
            # # ax.bar(flu.seg_lengths[current_seg_name]+0.01*flu.seg_lengths['seg%d_H1N1FULL' % 1],ylim,width=5,color='black')
            ax.tick_params(axis='both', which='major', labelsize=fontsize_xy_ticks_labels-3)
            ax.set_title( seg ,fontsize=fontsize_xy_label-3)
        return (fig,axs)






# influenza mutant class ...
class influenza_mutant(influenza):
    """docstring for influenza_mutant: 
    constructor takes an existing instance of the influenza object
    and modifies it in accordabce with the SNP requirements
    which can be muttion rate, genome positions, etc. (mutation rate - only now!)"""
    def __init__(self, wt, snp, fname='tmp.mutant.fasta'):
        # super(influenza_mutant, self).__init__()
        # self.arg = arg
        # make wt sequences mutable:
        local_genome_copy = dict([(seg,wt.genome[seg].tomutable())for seg in wt.genome])
        #
        # if snp is float, i.e. mutation rate ...
        #
        # let's store SNPs, so that we can access them afterwards ...
        self.snp = {}
        #
        for seg in wt.seglen:
            # number of SNPs to make
            snp_num = int(wt.seglen[seg]*float(snp))
            # unique indexes for SNP placement
            snp_idx = rnd.sample(xrange(wt.seglen[seg]),snp_num)
            #
            self.snp[seg] = np.asarray(snp_idx)
            # SNP placement ...
            for idx in snp_idx:
                local_genome_scopy[seg][idx] = rnd.choice(list( set('ATGC')-set(local_genome_copy[seg][idx]) ))
        ######################################################################################
        self.genome = dict([(seg,local_genome_copy[seg].toseq())for seg in local_genome_copy])
        #
        # # read influenza genome ...
        # self.genome = dict([(seq.id,seq.seq)for seq in SeqIO.parse(genome_fname,"fasta")])
        # # length of the segments ...
        self.seglen = copy.deepcopy(wt.seglen)
        # genome length ...
        self.genlen = sum(self.seglen.values())
        # number of segments ...
        self.segnum = len(self.genome)
        # # proteins ...
        # self.segcode = copy.deepcopy(wt.segcode)
        # read ORF data ...
        self.orf = copy.deepcopy(wt.orf)
        #
        #
        # ONE MORE THING!!!
        # mutated genome must be written to the disk! otherwise, we have no way of refolding it ...
        self.genome_fname = fname
        ####
        SeqIO.write([SeqRecord.SeqRecord(seq=self.genome[seg],id=seg,name='',description='') for seg in self.genome],self.genome_fname,'fasta')
        # now construction is complete, mutant genome is stored and can be accessed by RNA-folders; fname is a calss variable ...













