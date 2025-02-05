#!/usr/bin/python3

#NOTE! You need to setup a virtual environment in which to execute this script.
#See https://docs.google.com/document/d/1nn0T0OWEsQCBoCdaH6DSY69lQSbK3XnPlseyyQuU2Lc/edit for directions, but if you want to just copy-paste...
#   python3 -m venv ~/bin/merge_beagles_venv
#   source ~/bin/merge_beagles_venv/bin/activate
#   pip install dask
#	python -m pip install "dask[dataframe]" --upgrade

import argparse
from collections import OrderedDict
import dask.dataframe as dd

parser = argparse.ArgumentParser()
parser.add_argument('--beagles_list', '-b', help = 'A comma-separated list of gzipped beagle files (.beagle.gz). Do not include spaces.')
parser.add_argument('--missingness_cutoff', '-n', type = float, help = 'A float indicating the allowable site missingness; default = 0.25.', default = 0.25)
parser.add_argument('--maf_threshold', '-m', type = float, help = 'A float indicating the minimum minor allele frequency allowed; default = 0.05.', default = 0.05)
parser.add_argument('--output_filename', '-o', help = 'A filename (".beagle" will be added by the script); default = "merged".', default = "merged")
args = parser.parse_args()

beagles_list = args.beagles_list.split(",")

# dask to merge beagles
init_merged_beagle = dd.read_csv(beagles_list[0], sep='\t', compression='gzip', blocksize=None)
for b in beagles_list[1:]:
	next_df = dd.read_csv(b, sep='\t', compression='gzip', blocksize=None)
	init_merged_beagle = dd.merge(init_merged_beagle, next_df, how = 'outer', on = "marker")

# dask to write out (prepare for site-wise filtering)
init_merged_beagle.to_csv("init_merged.beagle", sep='\t', single_file=True, index = False)

# site-wise filtering (line-by-line processing)
## parse the header
with open("init_merged.beagle", 'r') as imb:
	init_merged_beagle_header_as_list = imb.readline().rstrip().split('\t')

num_inds = (len(init_merged_beagle_header_as_list) - (1 + len(beagles_list) * 2)) / 3
header_list = [] # a new header for the final (filtered) merged beagle
beagle_list_index = 0 # from left-to-right, the line will start with the first beagle the user passed
allele_cols_list = []
for i,c in enumerate(init_merged_beagle_header_as_list):
	if c == "marker":
		header_list = [c,"allele1", "allele2"]
	elif "allele1" in c:
		beagle_list_index += 1
		allele_cols_list.append(i)
	elif "allele2" in c:
		allele_cols_list.append(i)
	else:
		indiv_id = c.split('.')[0].split('_')[0] + "-" + str(beagle_list_index)
		header_list.append(indiv_id)
allele_cols_list.append(len(init_merged_beagle_header_as_list))

outfilename = args.output_filename + ".beagle"
with open(outfilename, 'w') as o:
	o.write('\t'.join(header_list) + '\n')
	with open("init_merged.beagle", 'r') as m:
		next(m)
		for raw_snp in m:
			#a switch to throw
			print_snp = True
			
			#variables
			snp = raw_snp.rstrip().split('\t')
#			print("SNP: " + ' '.join(snp))
			true_allele_orientation = []
			true_genotypes = []
			alleles_list = []
			for allele_index in allele_cols_list:
				if allele_index < len(snp):
					if snp[allele_index] != '':
						alleles_list.append(snp[allele_index])
			alleles_set = set(alleles_list)
#			print("Alleles set: " + ' '.join(alleles_set))
			if len(alleles_set) == 2:
				undirected_genotypes = OrderedDict()
				iterator = iter(alleles_list)
				allele_pairs = list(zip(iterator, iterator))
				for beagle_num, allele_pair in enumerate(allele_pairs):
					current_allele_index_index = beagle_num * 2 + 1
					current_allele_index = allele_cols_list[current_allele_index_index] + 1
					next_allele_index_index = beagle_num * 2 + 2
					next_allele_index = allele_cols_list[next_allele_index_index]
					snp_genotypes = snp[current_allele_index:next_allele_index]
					if '' not in snp_genotypes:
						allele_pair_tuple = tuple(allele_pair)
						if allele_pair in undirected_genotypes.keys():
							undirected_genotypes[allele_pair_tuple].extend(snp_genotypes)
						else:
							undirected_genotypes[allele_pair_tuple] = snp_genotypes
				#missingness filter (first draft)
				snp_string = ' '.join(snp)
#				print("snp string: " + snp_string)
				explicit_missing = snp_string.count("0.333333 0.333333 0.333333")
#				print("exp miss: " + str(explicit_missing))
				#implicit missing
				pos_and_alleles = ' '.join([' '.join(inner) for inner in allele_pairs]).count(".") + 1
#				print("pos and alleles count: " + str(pos_and_alleles))
				all_columns = snp_string.count(".")
				inds_with_data = (all_columns - pos_and_alleles) / 3
				implicit_missing = num_inds - inds_with_data
				perc_site_miss = (explicit_missing + implicit_missing) / num_inds
				#if snp passes missingness filter
				if perc_site_miss < args.missingness_cutoff:
					#check maf
					minor_allele_count = 0
					first_genotype = undirected_genotypes[list(undirected_genotypes.keys())[0]]
					first_hom_maj_index = first_hom_maj_count = first_hom_min_count = first_het_count = 0
					while first_hom_maj_index < len(first_genotype):
						first_het_index = first_hom_maj_index + 1
						first_hom_min_index = first_hom_maj_index + 2
						first_hom_maj_count += (float(first_genotype[first_hom_maj_index]) * 2)
						first_het_count += float(first_genotype[first_het_index])
						first_hom_min_count += (float(first_genotype[first_hom_min_index]) * 2)
						first_hom_maj_index += 3
					minor_allele_count = first_hom_min_count + first_het_count
					true_hom_maj_index = 0
					true_allele_orientation = list(undirected_genotypes.keys())[0]
					if len(undirected_genotypes.keys()) == 2:
						second_genotype = undirected_genotypes[list(undirected_genotypes.keys())[1]]
						second_hom_maj_index = second_hom_maj_count = second_hom_min_count = second_het_count = 0
						while second_hom_maj_index < len(second_genotype):
							second_het_index = second_hom_maj_index + 1
							second_hom_min_index = second_hom_maj_index + 2
							second_hom_maj_count += (float(second_genotype[second_hom_maj_index]) * 2)
							second_het_count += float(second_genotype[second_het_index])
							second_hom_min_count += (float(second_genotype[second_hom_min_index]) * 2)
							second_hom_maj_index += 3

						#check orientation of genotypes [0,1 vs 1,0]
						min_first_way = first_hom_min_count + second_hom_maj_count + first_het_count + second_het_count
						min_second_way = first_hom_maj_count + second_hom_min_count + first_het_count + second_het_count
						if min_second_way >= min_first_way:
							true_hom_maj_index = 0
							true_allele_orientation = list(undirected_genotypes.keys())[0]
							minor_allele_count = min_first_way
						else:
							true_hom_maj_index = 2
							true_allele_orientation = list(undirected_genotypes.keys())[1]
							minor_allele_count = min_second_way
					elif len(undirected_genotypes.keys()) > 2: #not biallelic (somehow):
						print_snp = False
						print("SNP not biallelic")

					#check maf	
					maf = minor_allele_count / (num_inds * 2)
					if maf >= args.maf_threshold:
						if true_hom_maj_index == 0:
							if len(undirected_genotypes.keys()) == 1:
								true_genotypes = first_genotype
							elif len(undirected_genotypes.keys()) == 2:
								while true_hom_maj_index < len(second_genotype):
									true_hom_min_index = true_hom_maj_index + 2
									#swap hom maj and hom min in second_genotypes
									second_genotype[true_hom_min_index], second_genotype[true_hom_maj_index] = second_genotype[true_hom_maj_index], second_genotype[true_hom_min_index]
									true_hom_maj_index += 3
								true_genotypes = first_genotype.extend(second_genotype)
						elif true_hom_maj_index == 2:
							true_hom_min_index = true_hom_maj_index - 2
					else: #if maf is too low
						print_snp = False
						print("SNP maf too low")
				else: #if missingness is too high
					print_snp = False
					print("SNP missingness too high")
			else: #not biallelic (somehow)
				print_snp = False
				print("SNP not biallelic")
			if print_snp == True:
				o.write('\t'.join([snp[0]] + list(true_allele_orientation) + true_genotypes) + '\n')