#!/usr/bin/python3
import argparse
import gzip
import os
from collections import OrderedDict

#Read in config file for the run
parser = argparse.ArgumentParser()
parser.add_argument('--ckpt_file', '-p', help = 'Please provide the checkpoint file created in step 0.')
parser.add_argument('--thinning_factor', '-t', type = float, help = 'Factor by which to thin the beagle file; that is, 1 of every t lines should be kept.', required = False)
args = parser.parse_args()

#Initialize run config variables with some default values
working_dir = None
scripts_dir = None
jobsout_dir = None
gls_dir = None
prefix = None
filtered_n = None
filtered_bamslist = None
ref_genome = None
wgph_beagle_files = None
email = None

#Parse the config file
with open(args.ckpt_file, 'r') as last_step_ckpt:
	for raw_ckpt_line in last_step_ckpt:
		ckpt_line = raw_ckpt_line.rstrip()
		ckpt_setting = ckpt_line.split('\t')
		if ckpt_setting[0] == "workingDIR":
			working_dir = ckpt_setting[1]
		elif ckpt_setting[0] == "scriptsDIR":
			scripts_dir = ckpt_setting[1]
		elif ckpt_setting[0] == "jobsoutDIR":
			jobsout_dir = ckpt_setting[1]
		elif ckpt_setting[0] == "glsDIR":
			gls_dir = ckpt_setting[1]
		elif ckpt_setting[0] == "usePREFIX":
			prefix = ckpt_setting[1]
		elif ckpt_setting[0] == "nIND-filtered":
			filtered_n = ckpt_setting[1]
		elif ckpt_setting[0] == "bamsLIST-filtered":
			filtered_bamslist_file = ckpt_setting[1]
		elif ckpt_setting[0] == "refgenomeFASTA":
			ref_genome = ckpt_setting[1]
		elif ckpt_setting[0] == "glsFILES-polymorphic-filtered":
			wgph_beagle_files = ckpt_setting[1].split(",")
		elif ckpt_setting[0] == "email":
			email = ckpt_setting[1]

#Make ld directory
ld_dir = working_dir + "linkage/"
if os.path.isdir(ld_dir) is not True:
	os.mkdir(ld_dir)

#Generate an array script (a job for every chromosome)
ld_array_input = scripts_dir + prefix + "_linkageARRAY_input.txt"
with open(ld_array_input, 'w') as ld:
	i = 0
	while i < len(wgph_beagle_files):
		ld.write(str(i + 1) + ":" + wgph_beagle_files[i] + "\n")
		i += 1

#Generate a script for ngsLD, which will calculate pairwise LD and prune linked sites
ld_script = scripts_dir + prefix + "_linkageARRAY.sh"
with open(ld_script, 'w') as s:
	s.write("#!/bin/bash\n\n")
	s.write("#SBATCH --cpus-per-task=10\n")
	s.write("#SBATCH --job-name=ngsLD\n")
	s.write("#SBATCH --output=" + jobsout_dir + prefix + "_ngsLD_%A-%a.out\n")
	s.write("#SBATCH --error=" + jobsout_dir + prefix + "_ngsLD_%A-%a.err\n")
	s.write("#SBATCH --time=7-00:00:00\n")
	s.write("#SBATCH --mail-type=FAIL\n")
	s.write("#SBATCH --mail-user=" + email + "\n")
	s.write("#SBATCH --array=1-" + str(len(wgph_beagle_files)) + "%8\n\n")

	s.write("module unload bio/ngsld/1.1.1 bio/angsd/0.933\n")
	s.write("module load bio/ngsld/1.1.1 bio/angsd/0.933\n\n")

	s.write("JOBS_FILE=" + ld_array_input + "\n")
	s.write("IDS=$(cat ${JOBS_FILE})\n\n")
	s.write("for sample_line in ${IDS}\n")
	s.write("do\n")
	s.write("""\tjob_index=$(echo ${sample_line} | awk -F ":" '{print $1}')\n""")
	s.write("""\tbeagle_file=$(echo ${sample_line} | awk -F ":" '{print $2}')\n\n""")
	s.write("\tif [[ ${SLURM_ARRAY_TASK_ID} == ${job_index} ]]; then\n")
	s.write("\t\tbreak\n")
	s.write("\tfi\n")
	s.write("done\n\n")

	s.write("base_filename=" + ld_dir + "$(echo $beagle_file | sed 's!^.*/!!')\n")

	if args.thinning_factor is not None:
		s.write("base_filename=${base_filename%.beagle.gz}_ngsLD-" + args.thinning_factor + "\n\n")
		#remove the first line and apply the thinning factor to the beagle file
		s.write("zcat ${beagle_file} | sed '1d' | awk 'NR % " + args.thinning_factor + " == 0' | gzip > " + \
		  "${base_filename}.beagle.gz\n\n")
		
		#get the associated pos file
		s.write("zcat ${base_filename}.beagle.gz | awk '{print $1}' | " + \
		  "sed -i 's/.1_/.1\t/g' > " + \
		  "${base_filename}.pos\n\n")
		
	else:
		s.write("base_filename=${base_filename%.beagle.gz}_ngsLD\n\n")
		#remove the first line of the beagle file
		s.write("zcat ${beagle_file} | sed '1d' | gzip > " + \
		  "${base_filename}.beagle.gz\n\n")
		
		#get the associated pos file
		s.write("zcat $base_filename}.beagle.gz | awk '{print $1}' | " + \
		  "sed -i 's/.1_/.1\t/g' > " + \
		  "${base_filename}.pos\n\n")

	s.write("zcat ${base_filename}.beagle.gz " + \
		"| cut -f 4- " + \
		"| gzip > ${base_filename}.beagle.gz\n\n")
	
	s.write('num_of_lines=$(< "${base_filename}.beagle.gz" zcat | wc -l)\n\n')

	s.write("ngsLD " + \
		"--geno ${base_filename}.beagle.gz " + \
		"--pos ${base_filename}.pos " + \
		"--probs " + \
		"--n_ind " + filtered_n + " " + \
		"--n_sites ${num_of_lines} " + \
		"--max_kb_dist 0 " + \
		"--n_threads 10 " + \
		"--out ${base_filename}.ld\n\n")
	
	s.write("prune_graph.pl " + \
		"--in_file ${base_filename}.ld " + \
		"--max_kb_dist 2 " + \
		"--min_weight 0.5 " + \
		"--out ${base_filename}-unlinked.id\n\n")

#Take the list of unlinked sites into ANGSD
#(not worried about doing a lot of filtering here because the beagle/maf that went in already represented filtered sites)
	# s.write("sed -i 's/:/\t/g' " + ld_dir + "${base_filename}-unlinked.id\n\n")
	# s.write("angsd sites index " + ld_dir + "${base_filename}-unlinked.id\n\n")
	# s.write("angsd -b " + filtered_bamslist + " " + \
	# 	"-ref " + ref_genome + " " + \
	# 	"-out " + gls_dir + "${base_filename}-unlinked " + \
	# 	"-nThreads 10 " + \
	# 	"-C 50 " + \
	# 	"-doCounts 1 " + \
	# 	"-GL 1 " + \
	# 	"-doGlf 2 " + \
	# 	"-doMajorMinor 3 " + \
	# 	"-doMaf 1 " + \
	# 	"-doPost 1 " + \
	# 	"-sites " + ld_dir + "${base_filename}-unlinked.id")

#Record new beagle and maf filenames for each chromosome
#unlinked_output = OrderedDict()
#for chrom in chrs_list:
#	unlinked_beagle = ld_dir + prefix + "_" + chrom + "_polymorphic-unlinked.beagle.gz"
#	unlinked_maf = ld_dir + prefix + "_" + chrom + "_polymorphic-unlinked.mafs.gz"
#	unlinked_output[chrom] = [unlinked_beagle, unlinked_maf]

#Update checkpoint file with the genotype likelihood file names
#with open(args.ckpt_file, 'a') as ckpt:
#	ckpt.write("all-chrs_mafsFILE-polymorphic\t" + mafs_file + "\n")
#	ckpt.write("linkageDIR\t" + ld_dir + "\n")
#	ckpt.write("ldFILE\t" + ld_dir + prefix + ".ld\n")
#	ckpt.write("unlinked-sitesFILE\t" + ld_dir + prefix + "-unlinked.id\n")
#	for chromosome, unlinked_files in unlinked_output.items():
#		ckpt.write("polymorphic-ldprunedFILES\t" + chromosome + "\t" + ",".join(unlinked_files) + "\n")