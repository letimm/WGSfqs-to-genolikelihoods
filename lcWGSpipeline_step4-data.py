#!/usr/bin/python3
import argparse
import gzip
import math
import os
import subprocess
from collections import OrderedDict
from math import sqrt

#Read in config file for the run
parser = argparse.ArgumentParser()
parser.add_argument('--ckpt_file', '-p', help = 'Please provide the checkpoint file created in step 0.')
parser.add_argument('--blacklist_inds', '-b', help = 'Please provide a file listing any individuals that should be removed from downstream analyses.')
args = parser.parse_args()

#Initialize run config variables with some default values
working_dir = None
scripts_dir = None
jobsout_dir = None
ref_genome = None
prefix = None
email = None
chrs_list = []
n_ind = None
endedness = None
blacklist_bams = []
full_bams_list = None
gls_filename = None

#Parse the config file to determine what's needed
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
		elif ckpt_setting[0] == "refgenomeFASTA":
			ref_genome = ckpt_setting[1]
		elif ckpt_setting[0] == "prefix":
			prefix = ckpt_setting[1]
		elif ckpt_setting[0] == "email":
			email = ckpt_setting[1]
		elif ckpt_setting[0] == "chrsLIST":
			chrs_list = ckpt_setting[1].split(",")
		elif ckpt_setting[0] == "nIND":
			n_ind = int(ckpt_setting[1])
		elif ckpt_setting[0] == "ENDEDNESS":
			endedness = ckpt_setting[1]
		elif ckpt_setting[0] == "bamsLIST-all":
			full_bams_list = ckpt_setting[1]

gls_dir = working_dir + "gls/"
if os.path.isdir(gls_dir) is not True:
	os.mkdir(gls_dir)

#Global angsd for gls
with open(args.blacklist_inds, 'r') as bb:
	for bb_id in bb:
		blacklist_bams.append(bb_id.rstrip())

filtered_bamslist_filename = working_dir + prefix + "_filtered_bamslist.txt"
with open(filtered_bamslist_filename, 'w') as fb:
	with open(full_bams_list, 'r') as b:
		for bam in b:
			blacklisted_status = False
			for black_bam in blacklist_bams:
				black_id = black_bam + "_"
				if black_id in bam:
					blacklisted_status = True
					n_ind -= 1
					break
			if blacklisted_status == False:
				fb.write(bam)

angsd_array_input = scripts_dir + prefix + "_angsdARRAY_input.txt"
with open(angsd_array_input, 'w') as i:
	iterator = 1
	for contig in chrs_list:
		i.write(str(iterator) + ":" + contig + "\n")
		iterator += 1

global_script = scripts_dir + prefix + "_globalARRAY.sh"
with open(global_script, 'w') as glb:
	glb.write("#!/bin/bash\n\n")
	glb.write("#SBATCH --cpus-per-task=10\n")
	glb.write("#SBATCH --time=0-20:00:00\n")
	glb.write("#SBATCH --job-name=global_" + prefix + "\n")
	glb.write("#SBATCH --output=" + jobsout_dir + prefix + "_global_%A-%a.out\n")
	glb.write("#SBATCH --mail-type=FAIL\n")
	glb.write("#SBATCH --mail-user=" + email + "\n")
	glb.write("#SBATCH --array=1-" + str(len(chrs_list)) + "%24\n\n")
	glb.write("module unload bio/angsd/0.933\n")
	glb.write("module load bio/angsd/0.933\n\n")

	glb.write("JOBS_FILE=" + angsd_array_input + "\n")
	glb.write("IDS=$(cat ${JOBS_FILE})\n\n")
	glb.write("for sample_line in ${IDS}\n")
	glb.write("do\n")
	glb.write("""\tjob_index=$(echo ${sample_line} | awk -F ":" '{print $1}')\n""")
	glb.write("""\tcontig=$(echo ${sample_line} | awk -F ":" '{print $2}')\n""")
	glb.write("\tif [[ ${SLURM_ARRAY_TASK_ID} == ${job_index} ]]; then\n")
	glb.write("\t\tbreak\n")
	glb.write("\tfi\n")
	glb.write("done\n\n")

	glb.write("angsd -b " + filtered_bamslist_filename + " " + \
		"-ref " + ref_genome + " " + \
		"-r ${contig}: " + \
		"-out " + gls_dir + prefix + "_${contig}_global " + \
		"-nThreads 10 " + \
		"-uniqueOnly 1 " + \
		"-remove_bads 1 " + \
		"-trim 0 " + \
		"-C 50 " + \
		"-minMapQ 15 " + \
		"-doCounts 1 " + \
		"-setminDepth " + str(n_ind) + " " + \
		"-setmaxDepth " + str(float(n_ind) * 20) + " " + \
		"-GL 1 " + \
		"-doGlf 2 " + \
		"-doMaf 1 " + \
		"-doMajorMinor 1 " + \
		"-doDepth 1 " + \
		"-dumpCounts 3")
	if endedness == "PE":
		glb.write(" -only_proper_pairs 1")

polymorphic_script = scripts_dir + prefix + "_polymorphicARRAY.sh"
with open(polymorphic_script, 'w') as plm:
	plm.write("#!/bin/bash\n\n")
	plm.write("#SBATCH --cpus-per-task=10\n")
	plm.write("#SBATCH --time=0-20:00:00\n")
	plm.write("#SBATCH --job-name=plm_" + prefix + "\n")
	plm.write("#SBATCH --output=" + jobsout_dir + prefix + "_polymorphic_%A-%a.out\n")
	plm.write("#SBATCH --mail-type=FAIL\n")
	plm.write("#SBATCH --mail-user=" + email + "\n")
	plm.write("#SBATCH --array=1-" + str(len(chrs_list)) + "%24\n\n")
	plm.write("module unload bio/angsd/0.933\n")
	plm.write("module load bio/angsd/0.933\n\n")

	plm.write("JOBS_FILE=" + angsd_array_input + "\n")
	plm.write("IDS=$(cat ${JOBS_FILE})\n\n")
	plm.write("for sample_line in ${IDS}\n")
	plm.write("do\n")
	plm.write("""\tjob_index=$(echo ${sample_line} | awk -F ":" '{print $1}')\n""")
	plm.write("""\tcontig=$(echo ${sample_line} | awk -F ":" '{print $2}')\n""")
	plm.write("\tif [[ ${SLURM_ARRAY_TASK_ID} == ${job_index} ]]; then\n")
	plm.write("\t\tbreak\n")
	plm.write("\tfi\n")
	plm.write("done\n\n")

	plm.write("angsd -b " + filtered_bamslist_filename + " " + \
		"-ref " + ref_genome + " " + \
		"-r ${contig}: " + \
		"-out " + gls_dir + prefix + "_${contig}_polymorphic " + \
		"-nThreads 10 " + \
		"-uniqueOnly 1 " + \
		"-remove_bads 1 " + \
		"-trim 0 " + \
		"-C 50 " + \
		"-minMapQ 15 " + \
		"-doCounts 1 " + \
		"-setminDepth " + str(n_ind) + " " + \
		"-setmaxDepth " + str(float(n_ind) * 20) + " " + \
		"-doGlf 2 " + \
		"-GL 1 " + \
		"-doMaf 1 " + \
		"-doMajorMinor 1 " + \
		"-minMaf 0.05 " + \
		"-SNP_pval 1e-10 " + \
		"-doDepth 1 " + \
		"-dumpCounts 3")
	if endedness == "PE":
		plm.write(" -only_proper_pairs 1")

#Update checkpoint file with the data file names
global_gls_files = []
polymorphic_gls_files = []
global_maf_files = []
polymorphic_maf_files = []
global_depths_files = []
polymorphic_depths_files = []
global_counts_files = []
polymorphic_counts_files = []

for chrom in chrs_list:
	global_gls_file = gls_dir + prefix + "_" + chrom + "_global.beagle.gz"
	global_gls_files.append(global_gls_file)
	polymorphic_gls_file = gls_dir + prefix + "_" + chrom + "_polymorphic.beagle.gz"
	polymorphic_gls_files.append(polymorphic_gls_file)
	global_maf_file = gls_dir + prefix + "_" + chrom + "_global.mafs.gz"
	global_maf_files.append(global_maf_file)
	polymorphic_maf_file = gls_dir + prefix + "_" + chrom + "_polymorphic.mafs.gz"
	polymorphic_maf_files.append(polymorphic_maf_file)
	global_depths_file = gls_dir + prefix + "_" + chrom + "_global.pos.gz"
	global_depths_files.append(global_depths_file)
	polymorphic_depths_file = gls_dir + prefix + "_" + chrom + "_polymorphic.pos.gz"
	polymorphic_depths_files.append(polymorphic_depths_file)
	global_counts_file = gls_dir + prefix + "_" + chrom + "_global.counts.gz"
	global_counts_files.append(global_counts_file)
	polymorphic_counts_file = gls_dir + prefix + "_" + chrom + "_polymorphic.counts.gz"
	polymorphic_counts_files.append(polymorphic_counts_file)

with open(args.ckpt_file, 'a') as ckpt:
	ckpt.write("glsDIR\t" + gls_dir + "\n")
	ckpt.write("blackLIST\t" + args.blacklist_inds + "\n")
	ckpt.write("nIND-filtered\t" + str(n_ind) + "\n")
	ckpt.write("bamsLIST-filtered\t" + filtered_bamslist_filename + "\n")
	ckpt.write("glsFILES-global\t" + ",".join(global_gls_files) + "\n")
	ckpt.write("glsFILES-polymorphic\t" + ",".join(polymorphic_gls_files) + "\n")
	ckpt.write("mafsFILES-global\t" + ",".join(global_maf_files) + "\n")
	ckpt.write("mafsFILES-polymorphic\t" + ",".join(polymorphic_maf_files) + "\n")
	ckpt.write("depthsFILES-global\t" + ",".join(global_depths_files) + "\n")
	ckpt.write("depthsFILES-polymorphic\t" + ",".join(polymorphic_depths_files) + "\n")
	ckpt.write("countsFILES-global\t" + ",".join(global_counts_files) + "\n")
	ckpt.write("countsFILES-polymorphic\t" + ",".join(polymorphic_counts_files) + "\n")

print("Step 4 has finished successfully! You will find two new scripts in ./scripts/:\n" + \
	global_script + " calculates genotype likelihoods across all sites on each chromosome (separately).\n" + \
	polymorphic_script + " calculates genotype likelihoods across all polymorphic sites on each chromosome (separately).\n")
print("Both scripts can run simultaneously.")
print("After they have run, you will have genotype likelihoods (gls) and allele frequencies (maf) for " + \
	"all sites in the genome (global) and putatively variable sites (polymorphic).")
