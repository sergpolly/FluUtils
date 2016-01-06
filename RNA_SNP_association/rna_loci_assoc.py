import sys
import re
import os
import string
import numpy as np
import pandas as pd
import subprocess as sub
# from Bio import SeqIO
import scipy.sparse as sparse
import scipy.stats as st
#
import flu_module as flu


#####################################################################################################
########### plotting ...
def plot_heatmap(data,pvals,shape,period,xticks,yticks,out_fname):
    # unpack some of the variables ...
    xdim,ydim = shape
    x_ticknames,x_top_ticknames = xticks
    y_ticknames = yticks
    # setup canves ...
    fig = plt.figure(figsize=(14,6))
    ax = fig.add_subplot(111)
    axt = ax.twiny()
    pic_extent = (0,xdim,0,ydim)
    # draw heatmap ...
    ax.imshow(data,interpolation='nearest',cmap = 'seismic',aspect='auto',vmin=-20,vmax=20,extent=pic_extent)
    # put p-values on heatmap ...
    for i in range(ydim):
        for j in range(xdim):
            ax.text(j+0.5,ydim-i-0.5,"%.f"%pvals[i,j],color='white',horizontalalignment='center',verticalalignment='center')
    # setup xticks & labels ...
    ax.set_xticks(np.arange(0.5,xdim+0.5,1))
    ax.set_yticks(np.arange(0.5,ydim+0.5,1))
    ax.set_xticklabels(x_ticknames)
    ax.set_yticklabels(y_ticknames)
    # xlims & xbounds ...
    ax.set_xlim((0,xdim))
    axt.set_xlim((0,xdim))
    ax.set_xbound(0,xdim)
    axt.set_xbound(0,xdim)
    # push ax_xlims to the axt ...
    ax_xlim = ax.get_xlim()
    axt.set_xlim(ax_xlim)
    #
    axt.set_xticks(np.arange(period-0.5,xdim,period))
    axt.set_xticklabels(x_top_ticknames)
    #
    ax.set_xlabel("loci group name")
    axt.set_xlabel(r"NP attraction, $\Delta g_{NP}$")
    ax.set_ylabel("H1N1 segment")
    #
    plt.savefig(out_fname)
#####################################################################################################


#####################################################################################################
def get_loci_structure_pattern(flu_obj,loci_of_interest,loci_title):
    # RANDOMIZATION BOOTSTRAP FUNCTION ...
    # THAT IS FOR P-VALUE CALCULATIONS ...
    # try randomization histograms ...
    def get_pval(value_observed,sample_size,profile,iterations=1000):
        # random sampling of the profile ...
        sample_profile  = lambda _profile: np.random.choice(_profile,size=sample_size,replace=False)
        # get fraction of structured features of the profile sample ...
        get_feature     = lambda _profile: get_struct(sample_profile(_profile))*100.0/sample_size
        # return percentile of the observed value (same as p-value) in randomly sampled (#iteration) distribution...
        return st.percentileofscore( [get_feature(profile) for _ in range(iterations)], value_observed )
    ######################################################
    get_bound   = lambda profile: (profile==0).nonzero()[0].size
    get_loops   = lambda profile: (profile==1).nonzero()[0].size
    get_paired  = lambda profile: (profile==2).nonzero()[0].size
    get_struct  = lambda profile: (profile!=0).nonzero()[0].size
    # store all expected vs observed fractions in a DataFrame ...
    bound = {"observed":[],"expected":[]}
    loops = {"observed":[],"expected":[]}
    paired = {"observed":[],"expected":[]}
    structured = {"observed":[],"expected":[]}
    long_loops = {"observed":[],"expected":[],"pval":[]}
    # for all segments extracts statisticall info on the loci - RNA structure association ...
    for seg in flu_obj.segments:
        # LOCI TO MAP ARE THE LOCATIONS OF SNP OF INTEREST ...
        loci_to_map = np.asarray(loci_of_interest.get_group(seg)['pos'])
        # get those profiles ...
        # just print out the number of mutations...
        print seg, " has %d loci of interest."%loci_to_map.size
        ####################################### ...
        pair_profile_seg = flu_obj.pair_profile[seg]
        long_loops_profile_seg = flu_obj.loops_profile[seg]
        ###################################################
        pair_profile_loci = pair_profile_seg[loci_to_map]
        long_loops_profile_loci = long_loops_profile_seg[loci_to_map]
        # bound,unbound,loops, etc. ratios observed for loci of interest ...
        bound_obs       = get_bound(pair_profile_loci)*100.0/pair_profile_loci.size
        loops_obs       = get_loops(pair_profile_loci)*100.0/pair_profile_loci.size
        paired_obs      = get_paired(pair_profile_loci)*100.0/pair_profile_loci.size
        struct_obs      = get_struct(pair_profile_loci)*100.0/pair_profile_loci.size
        long_loops_obs  = get_struct(long_loops_profile_loci)*100.0/long_loops_profile_loci.size
        # bound,unbound,loops, etc. ratios expected on average in seg ...
        bound_exp       = get_bound(pair_profile_seg)*100.0/flu_obj.seglen[seg]
        loops_exp       = get_loops(pair_profile_seg)*100.0/flu_obj.seglen[seg]
        paired_exp      = get_paired(pair_profile_seg)*100.0/flu_obj.seglen[seg]
        struct_exp      = get_struct(pair_profile_seg)*100.0/flu_obj.seglen[seg]
        long_loops_exp  = get_struct(long_loops_profile_seg)*100.0/flu_obj.seglen[seg]
        # fill in corresponding dictionary ...
        bound["observed"].append(bound_obs)
        loops["observed"].append(loops_obs)
        paired["observed"].append(paired_obs)
        structured["observed"].append(struct_obs)
        long_loops["observed"].append(long_loops_obs)
        #
        bound["expected"].append(bound_exp)
        loops["expected"].append(loops_exp)
        paired["expected"].append(paired_exp)
        structured["expected"].append(struct_exp)
        long_loops["expected"].append(long_loops_exp)
        # calculate long_loop structured feature p-value ...
        long_loops_pval = get_pval(long_loops_obs,
                                    long_loops_profile_loci.size,
                                    long_loops_profile_loci)
        # and append it ...
        long_loops["pval"].append(long_loops_pval)
    ############################
    #
    tmp = {}
    tmp['seg'] = sum([flu_obj.segments,flu_obj.segments],[])
    tmp['bound'] = sum([bound['observed'],bound['expected']],[])
    tmp['loops'] = sum([loops['observed'],loops['expected']],[])
    tmp['paired'] = sum([paired['observed'],paired['expected']],[])
    tmp['structured'] = sum([structured['observed'],structured['expected']],[])
    tmp['long_loops'] = sum([long_loops['observed'],long_loops['expected']],[])
    tmp['long_loops_pval'] = sum([long_loops['pval'],[None for _ in flu_obj.segments]],[])
    tmp['status'] = sum([['observed' for _ in flu_obj.segments],['expected' for _ in flu_obj.segments]],[])
    tmp['dgnp'] = sum([dg_NP for _ in flu_obj.segments],[dg_NP for _ in flu_obj.segments],[])
    tmp['loci'] = sum([loci_title for _ in flu_obj.segments],[loci_title for _ in flu_obj.segments],[])
    # store it to the big output ...
    return pd.DataFrame(tmp)






# MAIN PART STARTS HERE ...
#######################################################################################################

#NONVARIABLE THINGS GOES FIRST ...
#
# genome_fname = "../FluGenomes/H1N1FULL.fa"
genome_fname = "../FluGenomes/H1N1FULL.fa"
orf_fname = "../FluGenomes/xlator_data/H1N1FULL.orf"

# loci of interest must be in the form of DataFrame with the columns of different lengths ...
# or something else like that ...
loci_int = pd.read_csv('loci_of_interest.dat',index=False)
# columns = ['name','seg','pos']
loci_groups = loci_int.groupby('name')
loci_group_names = loci_groups.groups.keys()
loci_groups_num = len(loci_group_names)

########################################################
###                                                  ###
###      then we would vary folding conditions ...   ###
###                                                  ###
########################################################
# ['bound', 'dgnp', 'long_loops', 'loops_all', 'paired', 'structured', 'threshold_freq']
# # parameter to control RNA-NP non-specific attraction ... 
# dg_NP = float(sys.argv[1])
# threshold_freq = 1.0
enrichment_dat = pd.DataFrame({})
#############
dg_list = [-0.18,-0.2,-0.25,-0.3]
dg_list_num = len(dg_list)
###################
for dg_NP in dg_list:
    print "NP melting energy:",dg_NP
    ##################################
    # setup flu genome here ...
    h1n1 = flu.influenza(genome_fname,orf_fname)
    # calculate its RNA structure ...
    h1n1.calculate_rna_lMFE_optimized(dg_NP)
    ##################################
    for loci_group_name in loci_group_names:
        # print name ...
        print "loci of interest: ", loci_group_name
        loci = loci_groups.get_group(loci_group_name)
        loci_seg = loci[['seg','pos']].groupby('seg')
        #
        enrichment_tmp = get_loci_structure_pattern(h1n1, loci_seg, "%.3f"%threshold_freq)
        #
        enrichment_dat.append(enrichment_tmp,ignore_index=True)
####################################################
#  DEALING WITH THE ENRICHMENT DATA ...
####################################################
# group the enrichment data by 'observed' or 'expected' ...
enrichment_grouped = enrichment_dat.groupby("status")
# calculating the deviation of observed loci-rna association from expected average ...
# columns to get difference of ...
columns_delta = ['bound', 'long_loops', 'long_loops_pval', 'loops', 'paired', 'structured']
# difference ...
enrichment_observed = enrichment_grouped.get_group('observed').reset_index(drop=True)[columns_delta]
enrichment_expected = enrichment_grouped.get_group('expected').reset_index(drop=True)[columns_delta]
# difference calculation ...
enrichment_delta = enrichment_observed - enrichment_expected
# fill some supporting data ...
enrichment_delta['dgnp'] = np.asarray(enrichment_grouped.get_group('observed')['dgnp'])
enrichment_delta['threshold_freq'] = np.asarray(enrichment_grouped.get_group('observed')['threshold_freq'])
enrichment_delta['seg'] = np.asarray(enrichment_grouped.get_group('observed')['seg'])
####################################################################################


##################### PREPARATION FOR HEATMAP OUTPUT ###########################
heat_array_shape = (loci_groups_num*dg_list_num, h1n1.segnum)
# here one can change the association object:
# long_loops, or paired regions or simply structured parts ...
heat_data = enrichment_delta['long_loops'].reshape(heat_array_shape) # transpose before plotting ...
heat_pvals = enrichment_delta['long_loops_pval'].reshape(heat_array_shape) # transpose before plotting ...
##############################################################################
y_ticknames =enrichment_delta['seg'].reshape(heat_array_shape)[0]
x_ticknames = enrichment_delta['threshold_freq'].reshape(heat_array_shape)[:,0]
x_top_ticknames = enrichment_delta['dgnp'].reshape(heat_array_shape)[:,0]
x_top_ticknames = x_top_ticknames[range(dg_list_num-1,heat_array_shape[0],dg_list_num)]






















