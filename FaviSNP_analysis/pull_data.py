import sys
import re
import os
# import string
import numpy as np
# import pandas as pd
import subprocess as sub
from Bio import SeqIO
from Bio.Data import CodonTable
import pandas as pd
# import scipy.sparse as sparse
# import scipy.stats as st
# import flu_module as flu

from Bio import SeqIO
from Bio.Data import CodonTable


# Bib:
# /data/zeldovich/flustore/runs/RUN1234
# nucleoCountCumul.FULL
# /home/zeldovik/darpa/reassort-paper-FULL/flutraqGetBaseName.pl
# /home/zeldovik/darpa/reassort-paper-FULL /flutraqGetFORMOI.pl -r 1425


# cmd = ['RNALfold','-L',str(window)]
# RNALfold_process = sub.Popen(cmd,stdin=sub.PIPE,stdout=sub.PIPE)
# output,errput = RNALfold_process.communicate(local_seq_dict_fasta(seq_dict))
# local_struct_dict = local_parse_RNALfold(output)


# # run the optimizer ...
# cmd = "%s %s %s"%( os.path.join(bip_path,'run_bip'), fname_obj, fname_cst )
# result = sub.check_output(cmd,shell=True)
# # read last 'objective_size' lines from the result ...
# result_out = np.asarray( [ int(i) for i in result.strip().split('\n')[-objective_size:] ] )



# need this for some manipulations in .freq, .coverage and other flutraq data files ...
nt_index = dict( zip(list('ATGC'),range(4)) )

# genetic code ...
genetic_code = CodonTable.standard_dna_table.forward_table
stop_codons = dict([ (codon,'*') for codon in CodonTable.standard_dna_table.stop_codons ])
# make genetic code contain both 61 aacid codons and 3 stop codons ...
genetic_code.update(stop_codons)


# Kyte-Doolittle hydrophobicity scale ...
from Bio.SeqUtils import ProtParamData
KD = ProtParamData.kd
KD['*'] = 50. # just to make ->STOP, STOP-> to stand out... 


#########################################################################
# influenza class is needed mostly for ORF manipulations ...
class influenza(object):
    """influenza is a class describing all possible feature of a given influenza genome,
    including genome itself, lengths of segments, ORFs, RNA profiles etc."""
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





def flutraqGetBaseName(runid):
    # runid must be a single int value:
    runid = int(runid)
    #
    cmd = ["/home/zeldovik/darpa/reassort-paper-FULL/flutraqGetBaseName.pl","-r",str(runid)]
    output = sub.check_output(cmd)
    if output:
        return output.strip()
        # example: "RUN1234_H3N2_MDCK_pass002_barcodeTCCAAT\n"
    else:
        print "flutraqGetBaseName.pl returned nothing!"
        sys.exit(1)



def flutraqGetFORMOI(runid):
    # runid must be a single int value:
    runid = int(runid)
    #
    cmd = ["/home/zeldovik/darpa/reassort-paper-FULL/flutraqGetFORMOI.pl","-r",str(runid)]
    output = sub.check_output(cmd)
    if output:
        # return tuple(output.strip().split('\t')[1:])
        return pd.Series(output.strip().split('\t'),index=['runid','OS','Favi','Ribo','MOI'])
        # example "1234\tOS0\tF0\tR0\tM0.001\n"
    else:
        print "flutraqGetFORMOI.pl returned nothing!"
        sys.exit(1)



def flutraqGetDataName(runid,data_type='coverage'):
    # runid must be a single int value:
    runid = int(runid)
    #
    base_name = flutraqGetBaseName(runid)
    # base name example: "RUN1234_H3N2_MDCK_pass002_barcodeTCCAAT"
    strain = base_name.split('_')[1]+'FULL'
    # and the full base name is ...
    base_name_full = "%s.%s.blast.%s"%(base_name,strain,data_type)
    # the PATH "/data/zeldovich/flustore/runs/RUN1234/nucleoCountCumul.FULL"
    data_fname = "/data/zeldovich/flustore/runs/RUN%04d/nucleoCountCumul.FULL/%s"%(runid,base_name_full)
    return data_fname



def flutraqCheckDataIntegrity(runid_list):
    # get the data file name of the corresponding runid ...
    fname_list = [ flutraqGetDataName(runid) for runid in runid_list ]
    # check directories existence ...
    existence_list = [ os.path.exists(fname) for fname in fname_list ]
    # get existing runid-s only ...
    runid_exist = ','.join(runid for runid,existence in zip(runid_list,existence_list) if existence)
    # return quantity and runids themselves ...
    return sum(existence_list),runid_exist


# given a table with columns:
# ['seg','pos','ref','fA','fT','fG','fC','cover']
# we need several functions for statistics extraction ...

# Pandas(seg='seg1_H1N1FULL', pos=0, ref='A', fA=0.999791579825, fT=0.00020842017500000002, fG=0.0, fC=0.0, cover=4798)
def unroll_SNPs(snp_tuple):
    snp_record = []
    # unpack damn tuple ...
    seg,pos,ref,_,_,_,_,_ = snp_tuple
    # for all ATGC, but the reference one do:
    for nt in set(nt_index)-set(ref):
        # to access fA,fT,fG,fC we use numbered indexing of tuple:
        # index of fA: nt_index['A']+3, 3 because fA goes after seg,pos,ref !
        snp_record.append( (seg, pos, ref, "%s->%s"%(ref,nt), snp_tuple[nt_index[nt]+3]) )
    return snp_record

# using genomes and their ORFs infer codon substitutions ...
def getCodonSubs(snp_tuple,genomes):
    seg, pos, ref, snp, _, _, _, _, _, _ = snp_tuple
    # seg8_H1N1FULL
    _,strain = seg.split('_')
    the_genome = genomes[strain]
    orf = the_genome.orf[seg]
    segment = the_genome.genome[seg]
    # (orf_id, (codon_idx,codon_shift,codon_coord) )
    inorf = isin_orf(orf,pos)
    # codon substitutions ...
    codon_subs = [] 
    #
    # AT THIS POINT WE HAVE 2 SOURCES OF SEQUENCE INFORMATION: genome AND snp_tuple ITSELF:
    # make a simple check that reference nt matches.
    assert(segment[pos] == ref)
    #
    if not inorf:
        return None
    else:
        for product in inorf:
            aa_idx, codon_shift, codon_coord = inorf[product]
            # we use explicit coordinates for all 3 nucs in codon,
            # to handle within codon splice sites properly ...
            ref_codon = ''.join(  segment[_] for _ in codon_coord  )
            snp_codon = ''.join( (snp[-1] if (i==codon_shift) else segment[_]) for i,_ in enumerate(codon_coord) )
            codon_subs.append("%s->%s"%(ref_codon,snp_codon))
        return ','.join(codon_subs)

def translateCStoAAS(codon_snps):
    # 'ATG->GTC,CAT->AAA' translated to 'M->V,H->K'
    return ','.join('->'.join(genetic_code[c] for c in  cs.split('->')) for cs in codon_snps.split(','))

def getdKDforAAS(aa_subs):
    # 'M->V,H->K' to 'KD(M)-KD(V),KD(H)-KD(K)'
    get_dKD = lambda aa1,aa2: KD[aa2]-KD[aa1]
    return ','.join( "%.3f"%get_dKD( *aas.split('->') ) for aas in aa_subs.split(','))


#######################################################################
# load H1N1 and H3N2 genomes for ORF and amino acid level analysis ...
genomes = {}
# genome_fname = "../FluGenomes/H1N1FULL.fa"
genome_fname = "../FluGenomes/H1N1FULL.fa"
orf_fname = "../FluGenomes/xlator_data/H1N1FULL.orf"
#
h1n1 = influenza(genome_fname,orf_fname)
genomes['H1N1FULL'] = h1n1

# genome_fname = "../FluGenomes/H3N2FULL.fa"
genome_fname = "../FluGenomes/H3N2FULL.fa"
orf_fname = "../FluGenomes/xlator_data/H3N2FULL.orf"
#
h3n2 = influenza(genome_fname,orf_fname)
genomes['H3N2FULL'] = h3n2
#######################################################################



# upload data ...
run_data = pd.read_csv('runid.csv')

# because there are multiple instances of 'runid' for some entries,
# let's choose workable ones first!
run_data['num_runid_exist'],run_data['runid_exist'] = zip(*run_data['runid'].str.split('/').apply(flutraqCheckDataIntegrity))



# if we have redundancy, let human decide what to do ...
if (run_data['num_runid_exist']!=1).any():
    print "There are redundant runids! Impossible to proceede ..."
    sys.exit(1)
else:
    print "Data seem well-defined! Do analysis!"


# GET INFO FOR EVERY INSTANCE ...
# flutraqGetFORMOI(runid)
# example "1234\tOS0\tF0\tR0\tM0.001\n"
run_info = run_data['runid_exist'].apply(flutraqGetFORMOI)
# merge data tables ...
run_data = run_data.merge(run_info,left_on='runid_exist',right_on='runid',suffixes=['_all','_duplicate'])
#
# fix Ribo data, 'R' is apparently no info thing, should be changed to 'R0', implying no Ribo added ...
run_data['Ribo'] = run_data['Ribo'].replace('R','R0')



# for each runid, corresponding data files 'coverage' or 'count-all' or whatever, weights 
# less than or about 1MB each, so we can safely and easily upload them into memory ...
# it's ~100*(1Mb+1Mb) ~200Mb, both for coverage and for frequencies ...
# we've expanded that significantly since then: ~5X or something ...
snp_info = {}
for runid in run_data['runid_exist']:
    print "Performing analysis for %s ..."%str(runid) 
    # coverage uploading ...
    coverage_fname = flutraqGetDataName(runid,data_type='coverage')
    tmp_coverage = pd.read_csv(coverage_fname,sep=' ',names=['seg','pos','ref','cover'])
    # frequencies uploading ...
    freq_fname = flutraqGetDataName(runid,data_type='ntfreq')
    tmp_freq = pd.read_csv(freq_fname,sep=' ',names=['seg','pos','ref','fA','fT','fG','fC'])
    # merge these 2 tables ...
    tmp_snp_info = tmp_freq.merge(tmp_coverage,on=['seg','pos','ref'])
    # unroll SNP information for convenience ...
    unrolled_SNP = sum(( unroll_SNPs(row) for row in tmp_snp_info.itertuples(index=False) ),[])
    tmp_unrolled_SNP_df = pd.DataFrame(unrolled_SNP,columns=['seg','pos','ref','snp','freq'])
    # merge and store ...
    tmp_snp_info = tmp_unrolled_SNP_df.merge(tmp_snp_info,on=['seg','pos','ref'])
    #
    # # WE SHOULD INFER AND WRITE DOWN CORRESPONDING AMINO ACID SUBSTITUTION HERE AS WELL ...
    # # assuming genomes are loaded and we know the strain here ...
    tmp_snp_info['snp_codon'] = [ getCodonSubs(row,genomes=genomes) for row in tmp_snp_info.itertuples(index=False) ]
    #
    # Then we'd need to get AA substitutions based on codons and correspondingly the dKD for those ...
    # TAKES ~15 MINUTES ...
    # ...
    snp_info[runid] = tmp_snp_info



for runid in run_data['runid_exist']:
    snp_info[runid].to_csv("%s.csv"%str(runid),index=False)





# # runid = '1359'
# unrolled_SNP = sum(( unroll_SNPs(row) for row in snp_info[runid].itertuples(index=False) ),[])








