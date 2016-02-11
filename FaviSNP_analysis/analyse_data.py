import sys
import re
import os
# import string
import numpy as np
import subprocess as sub
from Bio import SeqIO
from Bio.Data import CodonTable
import pandas as pd
# import scipy.sparse as sparse
# import scipy.stats as st
# import flu_module as flu
from Bio import SeqIO
from Bio.Data import CodonTable
import itertools
import matplotlib.pyplot as plt



############################################################################################


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

############################################################################################


favi_data_path = "~/DATA/FaviProcessedData"


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


def translateCStoAAS(codon_snps):
    # 'ATG->GTC,CAT->AAA' translated to 'M->V,H->K'
    return ','.join('->'.join(genetic_code[c] for c in  cs.split('->')) for cs in codon_snps.split(','))

def getdKDforAAS(aa_subs):
    # 'M->V,H->K' to 'KD(M)-KD(V),KD(H)-KD(K)'
    get_dKD = lambda aa1,aa2: KD[aa2]-KD[aa1]
    return ','.join( "%.3f"%get_dKD( *aas.split('->') ) for aas in aa_subs.split(','))


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
# snp_info[runid] = tmp_snp_info


for runid in run_data['runid_exist']:
    snp_info[runid] = pd.read_csv( os.path.join(favi_data_path,"%s.csv"%str(runid)) )




def getSNPcount(df,fst,dst):
    return df[(df['freq']>fst)&(df['cover']>dst)].shape[0]


def getSNPdensity(df,fst,dst):
    count = df[(df['freq']>fst)&(df['cover']>dst)].shape[0]
    count_tot = df[(df['cover']>dst)].shape[0]
    return count/float(count_tot)


def get_drug_type(drug_tuple):
    ost,favi,ribo = drug_tuple
    ost = float( ost.lstrip('OS') )
    favi = float( favi.lstrip('F') )
    ribo = float( ribo.lstrip('R') )
    if (ost<1e-4)and(favi<1e-4)and(ribo<1e-4):
        return "NOdrug"
    elif ost>1e-4:
        return "OSelt"
    elif favi>1e-4:
        return "Favi"
    elif ribo>1e-4:
        return "Ribo"
    else:
        print "Supposedly, there are no mixtures of drugs in the data ..."
        print "Termination!"
        sys.exit(1)

#####################################
# # group run_data by drug presence/concentration 'OS' 'Favi' 'Ribo'
# run_drug_group = run_data.groupby(by=['OS','Favi','Ribo'])
# runid_grouped = run_drug_group['runid_exist']
run_data['drug'] = run_data[['OS','Favi','Ribo']].apply(get_drug_type,axis=1)
run_drug_group = run_data.groupby(by='drug')
runid_grouped = run_drug_group['runid_exist']
colors = itertools.cycle(["r", "b", "g"])


def bar_grouped_drug(drug_group,fst,dst):
    # for rects1 = ax.bar(ind, menMeans, width, color='r', yerr=menStd)
    x_ind = 0
    bar_width = 0.3
    x_ticks = []
    for drug, runid_exist in drug_group:
        tick_shift = 0
        color = next(colors)
        for runid in runid_exist:
            bar_height = getSNPcount(snp_info[runid],fst,dst)
            plt.bar(x_ind,bar_height,bar_width,color=color)
            x_ind += bar_width
            tick_shift += bar_width
        # collect x tciks just in case ...
        x_ticks.append( (drug,x_ind-0.5*tick_shift) )
        x_ind += 3.0*bar_width # make a gap between the groups ...
    #
    labels,ticks_pos = zip(*x_ticks)
    ax = plt.gca()
    ax.set_xticks(ticks_pos)
    ax.set_xticklabels(labels)
    ax.set_ylabel('SNP count')
    ax.set_xlabel('drug group')
    ax.set_title('made for d>%.1f and f>%.4f'%(dst,fst))
    plt.show()



fst_list = [0.005,0.01,0.3]
dst_list = [100,200,500,1000]
# let's create a pandas PANEL  of SNP counts with items - runids, major axis d_star and minor axis f_star ...
to_panel = {}
for runid in run_data['runid_exist']:
    to_panel[runid] = {}
    for dst in dst_list:
        to_panel[runid][dst] = {}
        for fst in fst_list:
            value = getSNPcount(snp_info[runid],fst,dst)
            # print runid,fst,dst,value
            # to_panel[runid][dst].append( value )
            to_panel[runid][dst][fst] = value

# make it a panel with corresponding indexes ...
rdf_panel = pd.Panel(to_panel)


# rdf_panel['981'].loc[0.010,100]







