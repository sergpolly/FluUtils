import sys
import re
import os
# import string
import numpy as np
# import pandas as pd
# import subprocess as sub
from Bio import SeqIO
from Bio.Data import CodonTable
import matplotlib.pyplot as plt

# Kyte-Doolittle hydrophobicity scale ...
from Bio.SeqUtils import ProtParamData
KD = ProtParamData.kd

# 12 SNP types ...
snp_types = ['->'.join([n1,n2]) for n1 in "ATGC" for n2 in "ATGC" if n1!=n2]

# initialize SNP storage ...
snp = {}
for snp_type in snp_types:
	snp[snp_type] = []


# seg1_H1N1FULL_A1T
def get_dna_snp_type(snp):
	snp = snp.strip().split('_')[-1]
	dna_letters = '|'.join('ATGC')
	snp_match = re.match("(%s)\d+(%s)"%(dna_letters,dna_letters),snp)
	if not snp_match:
		raise ValueError("dna SNP code is invalid %s ..."%snp)
	else:
		n_from	= snp_match.group(1)
		n_to	= snp_match.group(2)
	return '->'.join([n_from,n_to])


def get_aa_subs_pair(snp):
	snp = snp.strip().split('_')[-1]
	aa_letters = '|'.join(KD.keys())
	snp_match = re.match("(%s)\d+(%s)"%(aa_letters,aa_letters),snp)
	if not snp_match:
		raise ValueError("aa SNP code is invalid %s ..."%snp)
	else:
		aa_from	= snp_match.group(1)
		aa_to	= snp_match.group(2)
	return (aa_from,aa_to)


def get_dKD(aa_from,aa_to):
	return KD[aa_to]-KD[aa_from]



def plot_dKD(snp_type,all_dKD,bins=50):
	data = np.asarray(all_dKD[snp_type])
	plt.clf()
	plt.hist(data,bins=bins)
	ax = plt.gca()
	the_mean = data.mean()
	the_std = data.std()
	ax.text(0.75,0.9,"mean: %.3f\nstd:   %.3f"%(the_mean,the_std),transform=ax.transAxes)
	ax.set_title(snp_type)
	ax.set_xlabel('delta_KD')
	plt.savefig(snp_type+'.pdf')



fname = "snp_trans.cache"
with open(fname,'r') as fp:
	# accumulate interesting SNP ...
	for line in fp.readlines():
		items = line.strip().split()
		# unpack data ...
		dna_snp, aa_snps = items[0], items[1:]
		# filter NORF-s
		if aa_snps == ['NORF']:
			pass
		else:
			# dna snp type ...
			snp_type = get_dna_snp_type(dna_snp)
			for aa_snp in aa_snps:
				# filter aa mutations involving STOP(*) codons
				if aa_snp.find('*' )== -1:
					delta_KD = get_dKD(*get_aa_subs_pair(aa_snp))
					snp[snp_type].append(delta_KD)


###### now snp dictionary holds all amino acid substitutions (STOP codons excluded) grouped by snp_type ...
# snp holds dKDs grouped by dna snp types ...





for snp_type in snp_types:
	plot_dKD(snp_type,snp,bins=30)






