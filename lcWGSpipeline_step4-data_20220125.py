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
gls_dir = None
prefix = None
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
		elif ckpt_setting[0] == "glsDIR":
			gls_dir = ckpt_setting[1]
		elif ckpt_setting[0] == "prefix":
			prefix = ckpt_setting[1]
		elif ckpt_setting[0] == "chrsLIST":
			chrs_list = ckpt_setting[1].split(",")
		elif ckpt_setting[0] == "nIND":
			n_ind = int(ckpt_setting[1])
		elif ckpt_setting[0] == "ENDEDNESS":
			endedness = ckpt_setting[1]
		elif ckpt_setting[0] == "bamsLIST-all":
			full_bams_list = ckpt_setting[1]

#Global angsd for gls
with open(args.blacklist_inds, 'r') as bb:
	for bb_id in bb:
		blacklist_bams.append(bb_id.rstrip())

filtered_bamslist_filename = working_dir + "filtered_bamslist.txt"
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

gls_script = scripts_dir + prefix + "_glsARRAY.sh"
with open(gls_script, 'w') as gls:
	gls.write("#!/bin/bash\n\n")
	gls.write("#SBATCH --cpus-per-task=10\n")
	gls.write("#SBATCH --time=7-00:00:00\n")
	gls.write("#SBATCH --job-name=gls_" + prefix + "\n")
	gls.write("#SBATCH --output=" + jobsout_dir + prefix + "_gls_%A-%a.out\n")
	gls.write("#SBATCH --array=1-" + str(len(chrs_list)) + "%24\n\n")
	gls.write("module unload bio/angsd/0.933\n")
	gls.write("module load bio/angsd/0.933\n\n")
	gls.write("JOBS_FILE=" + angsd_array_input + "\n")
	gls.write("IDS=$(cat ${JOBS_FILE})\n\n")
	gls.write("for sample_line in ${IDS}\n")
	gls.write("do\n")
	gls.write("""\tjob_index=$(echo ${sample_line} | awk -F ":" '{print $1}')\n""")
	gls.write("""\tcontig=$(echo ${sample_line} | awk -F ":" '{print $2}')\n""")
	gls.write("\tif [[ ${SLURM_ARRAY_TASK_ID} == ${job_index} ]]; then\n")
	gls.write("\t\tbreak\n")
	gls.write("\tfi\n")
	gls.write("done\n\n")
	gls.write("angsd -b " + filtered_bamslist_filename + " " + \
		"-ref " + ref_genome + " " + \
		"-out " + gls_dir + prefix + "-${contig}_global " + \
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
		"-doMajorMinor 1 " + \
		"-doDepth 1 " + \
		"-dumpCounts 3")
	if endedness == "PE":
		gls.write(" -only_proper_pairs 1")

maf_script = scripts_dir + prefix + "_mafsARRAY.sh"
with open(maf_script, 'w') as maf:
	maf.write("#!/bin/bash\n\n")
	maf.write("#SBATCH --cpus-per-task=10\n")
	maf.write("#SBATCH --time=7-00:00:00\n")
	maf.write("#SBATCH --job-name=maf_" + prefix + "\n")
	maf.write("#SBATCH --output=" + jobsout_dir + prefix + "_maf_%A-%a.out\n")
	maf.write("#SBATCH --array=1-" + str(len(chrs_list)) + "%24\n\n")
	maf.write("module unload bio/angsd/0.933\n")
	maf.write("module load bio/angsd/0.933\n\n")
	maf.write("JOBS_FILE=" + angsd_array_input + "\n")
	maf.write("IDS=$(cat ${JOBS_FILE})\n\n")
	maf.write("for sample_line in ${IDS}\n")
	maf.write("do\n")
	maf.write("""\tjob_index=$(echo ${sample_line} | awk -F ":" '{print $1}')\n""")
	maf.write("""\tcontig=$(echo ${sample_line} | awk -F ":" '{print $2}')\n""")
	maf.write("\tif [[ ${SLURM_ARRAY_TASK_ID} == ${job_index} ]]; then\n")
	maf.write("\t\tbreak\n")
	maf.write("\tfi\n")
	maf.write("done\n\n")
	maf.write("angsd -b " + filtered_bamslist_filename + " " + \
		"-ref " + ref_genome + " " + \
		"-out " + gls_dir + prefix + "-${contig}_global " + \
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
		"-doMaf 1 " + \
		"-doMajorMinor 1 " + \
		"-minMaf 0.05 " + \
		"-SNP_pval 1e-10 " + \
		"-doDepth 1 " + \
		"-dumpCounts 4")
	if endedness == "PE":
		maf.write(" -only_proper_pairs 1")
	maf.write("\n\n")
	maf.write("gunzip -c " + gls_dir + prefix + "-${contig}_global " + \
		"| cut -f 1,2,3,4 " + \
		"| tail -n +2 > " + \
		gls_dir + prefix + "-${contig}_SNPs.txt\n\n")
	maf.write("angsd sites index " + gls_dir + prefix + "-${contigs}_SNPs.txt\n\n")
	maf.write("angsd -b " + filtered_bamslist_filename + " " + \
		"-ref " + ref_genome + " " + \
		"-out " + gls_dir + prefix + "-${contig}_polymorphic " + \
		"-nThreads 10 " + \
		"-C 50 " + \
		"-doCounts 1 " + \
		"-GL 1 " + \
		"-doGlf 2 " + \
		"-doMajorMinor 3 " + \
		"-doMaf 1 " + \
		"-doPost 1 " + \
		"-sites " + gls_dir + prefix + "-${contig}_SNPs.txt")

#Update checkpoint file with the data file names
gls_global_files = []
gls_polymorphic_files = []
maf_global_files = []
maf_polymorphic_files = []
depths_files = []
counts_files = []

for chrom in chrs_list:
	gls_global_file = gls_dir + prefix + "-" + chrom + "_global.beagle.gz"
	gls_global_files.append(gls_global_file)
	gls_polymorphic_file = gls_dir + prefix + "-" + chrom + "_polymorphic.beagle.gz"
	gls_polymorphic_files.append(gls_polymorphic_file)
	maf_global_file = gls_dir + prefix + "-" + chrom + "_global.mafs.gz"
	maf_global_files.append(maf_global_file)
	maf_polymorphic_file = gls_dir + prefix + "-" + chrom + "_polymorphic.mafs.gz"
	maf_polymorphic_files.append(maf_polymorphic_file)
	depths_file = gls_dir + prefix + "-" + chrom + "_global.pos.gz"
	depths_files.append(depths_file)
	counts_file = gls_dir + prefix + "-" + chrom + "_global.counts.gz"
	counts_files.append(counts_file)

with open(args.ckpt_file, 'a') as ckpt:
	ckpt.write("blackLIST\t" + args.blacklist_inds + "\n")
	ckpt.write("nIND-filtered\t" + str(n_ind) + "\n")
	ckpt.write("bamsLIST-filtered\t" + filtered_bamslist_filename + "\n")
	ckpt.write("glsFILES-global\t" + ",".join(gls_global_files) + "\n")
	ckpt.write("glsFILES-polymorphic\t" + ",".join(gls_polymorphic_files) + "\n")
	ckpt.write("mafsFILES-global\t" + ",".join(maf_global_files) + "\n")
	ckpt.write("mafsFILES-polymorphic\t" + ",".join(maf_polymorphic_files) + "\n")
	ckpt.write("depthsFILES\t" + ",".join(depths_files) + "\n")
	ckpt.write("countsFILES\t" + ",".join(counts_files) + "\n")

print("Step 4 has finished successfully! You will find two new scripts in ./scripts/: " + \
	gls_script + " and " + maf_script + ".")
print("These scripts can be run simultaneously.")
print("After they have run, you will have genotype likelihoods (gls) and allele frequencies (maf) for " + \
	"all sites in the genome (global) and putatively variable sites (polymorphic).")
print("Congratulations! This represents the end of data assembly.")
