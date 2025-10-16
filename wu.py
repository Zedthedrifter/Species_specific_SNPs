#!/usr/local/bin/python
# 2021.04.21
# Wu Huang<whuang@rbge.org.uk>
# Input: vcf file
# Output: AIL considering hetero info (multiple alleles suited)

from __future__ import division
import sys
import gzip
import copy
import argparse

## files
parser = argparse.ArgumentParser(description='\
extract ancestral informative SNPs from vcf files for target capture data\
required arguments:\
[-f | -d | -n | -s | -o]')
parser.add_argument('-v','--vcf', help='vcf file with GT and DP info')
parser.add_argument('-n','--names', help='sample ID to species names corresponding file. Format - ID scientific_name')
parser.add_argument('-sp','--selectedspps', help='a list of spps names that are multiple sampled (scientific_name)')
parser.add_argument('-sm','--selectedsample', help='a list of sample names that belongs to target genus (ID)')
parser.add_argument('-o','--output',help='prefix of output files, can include the path to the directory for the output files')
args = parser.parse_args()

all_samples_to_spps_dict = {}

samples_in_selected_species = {}
count_bp = {}
for i in open(args.selectedspps):
	print(i)
	i = i.strip().split()
	samples_in_selected_species[i[0]] = []
	count_bp[i[0]] = 0

selected_samples = []
for i in open(args.selectedsample):
	print(i)
	i = i.strip().split()
	selected_samples.append(i[0])

for i in open(args.names):
	i = i.strip().split()
	print(i)
	all_samples_to_spps_dict[i[0]] = i[1]
	if i[1] in samples_in_selected_species:
		samples_in_selected_species[i[1]].append(i[0])

#initiate the two output files
AIL_out_file = open(args.output+'.AIL.list','w')
AIL_out_file.write('Spp\tChr\tPos\tType\tSpp_A_Allele\tRest_ind_allele\tA_max_AF\tAF_A_list\tAF_ALL_list\trest_selected_species_AF_list\n')
COUNT = open(args.output+'.count','w')
COUNT.write('Spp\tValid_Basepairs\n')

#calculate allele freq
def allele_freq(species_alleles):
	freq = [0,0,0,0]
	for ind_allele in species_alleles:
		if '.' in ind_allele:
			continue
		else:
			ind_allele = ind_allele.split('/')
			for allele in ind_allele:
				freq[int(allele)] += 1
	return [freq,sum(freq)]

# for every SNP locus, put all individuals's allele in a list: all_alleles = ['','',0/0,1/0,'',1/0,0/0,'',0/1,1/0,'0/2','','',...]
if args.vcf.endswith(".gz"):
	opener = gzip.open
else:
	opener = open

for i in opener(args.vcf, "rt"):
	if i.startswith('#CHROM'):
		position_in_vcf = i.strip().split()
		Nsample=len(position_in_vcf) #this is 9+ sample number? 
		print(Nsample)
## position_in_vcf = [#CHROM, POS, ID,...]
	else:
		if not i.startswith('#'): #each line in vcf, i.e. each SNP
			all_alleles = ['']*Nsample
			i = i.strip().split()
			for order in range(9,Nsample): #starting from 9th column
				flags = i[order].split(':')
				GT = flags[0] #genotype
				all_alleles[order] = GT
			#print(all_alleles) #genotype of all samples
## for every species, put its individuals into a list, the rest another list, create a list for each rest species and put them into a dictionary-{rest_species_allele_dict}
			for A in samples_in_selected_species:
				rest_alleles = []
				species_alleles = []
				rest_individual_list, rest_species_name, rest_species_allele_dict = selected_samples[:], copy.deepcopy(samples_in_selected_species),{}
				del rest_species_name[A]
				for individuals in samples_in_selected_species[A]:
					rest_individual_list.remove(individuals)
					#print(position_in_vcf.index(individuals))
					species_alleles.append(all_alleles[position_in_vcf.index(individuals)])
				print(species_alleles)                                                                                 
## species_alleles ['1/1', '1/1', '1/1', '1/1']
				for rest_individuals in rest_individual_list:
					rest_alleles.append(all_alleles[position_in_vcf.index(rest_individuals)])
				freq_A, freq_ALL = allele_freq(species_alleles), allele_freq(rest_alleles)
				rest_selected_species_AF_list = [[]for a in range(4)]
				for B in rest_species_name:
					B_species_alleles = []
					for B_individuals in rest_species_name[B]:
						B_species_alleles.append(all_alleles[position_in_vcf.index(B_individuals)])
					rest_species_allele_dict[B] = B_species_alleles
					freq_B = allele_freq(B_species_alleles)
					if freq_B[1] != 0:
						for n in range(0,4):
							rest_selected_species_AF_list[n].append(freq_B[0][n]/freq_B[1])
					else:
						continue
# judge if alleles_frequency in species A is significantly different from all the rest species respetively and in all.
## retaining loci by thresholds
				if freq_A[1] >= 2: # for this locus, species A must have at least 2 inds have data
					count_bp[A] += 1
					if len(rest_selected_species_AF_list[0]) >= 2: ## for this locus, alleles in at least 2 other groups are not missing
						AF_A_list = []
						AF_ALL_list =[]
						for n in range(0,4):
							AF_A_list.append(freq_A[0][n]/freq_A[1])
							AF_ALL_list.append(freq_ALL[0][n]/freq_ALL[1])
						max_index=AF_A_list.index(max(AF_A_list))
# stage 1: If AF in this group is higher than a specific threshold (say 87.5%), progresses to stage 2 ->
						if max(AF_A_list) >= 0.875:
# stage 2: if AF in all the rest groups(aggregated) is lower than a threshold (10%), progress to stage 3 ->
							if AF_ALL_list[max_index] <= 0.1:
								print(i[0],i[1])
								print('A',A,species_alleles)
								print(freq_A)
								print('ALL',rest_alleles)
								print('B',rest_species_allele_dict,rest_selected_species_AF_list)
# stage 3: if AFs in any other groups (separated) are lower than a threshold (12.5%, at most one allele amongst 8 alleles of 4 individuals), take this locus as valid one with allele frequency difference
								if max(rest_selected_species_AF_list[max_index]) <= 0.125:
									if max(AF_A_list) == 1.0 and AF_ALL_list[max_index] == 0.0: # IF fixed, TAG the SNP with 'SSA'
										AIL_out_file.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n'.format(A,i[0], i[1],'ssa',freq_A[0], freq_ALL[0], round(max(AF_A_list),2), AF_A_list,AF_ALL_list, rest_selected_species_AF_list))
									else:
										AIL_out_file.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n'.format(A, i[0], i[1],'.', freq_A[0], freq_ALL[0], round(max(AF_A_list),2), AF_A_list,AF_ALL_list, rest_selected_species_AF_list))
								else:
									continue
							else:
								continue
for A in count_bp:
	COUNT.write('{0}\t{1}\n'.format(A,count_bp[A]))

AIL_out_file.close()
