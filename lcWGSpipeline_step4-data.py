#!/usr/bin/python3
import argparse
import math
import os
import shutil
from collections import OrderedDict
from math import sqrt

#Read in config file for the run
parser = argparse.ArgumentParser()
parser.add_argument('--ckpt_file', '-p', help = 'Please provide the checkpoint file created in step 0.')
parser.add_argument('--blacklist_inds', '-b', help = 'Please provide a file listing any individuals that should be removed from downstream analyses.')
parser.add_argument('--global_gls', action = 'store_true', help = 'Set this flag if you would like to calculate genotype likelihoods for all sites across the genome (if this flag is not set, gls will only be calculated for SNPs).')
parser.add_argument('--quality_val', '-q', default = '15', help = 'Please provide a value for minQ and minMapQ. Default = 15')
parser.add_argument('--minInd', '-i', help = 'You can use this flag to specify a minInd parameter.', required = False)
parser.add_argument('--minDepth_factor', '-d', type = float, default = '1', help = 'Please provide a factor to define minDepth: minDepth_factor * N = minDepth. Default = 1')
parser.add_argument('--maxDepth_factor', '-D', type = float, default = '5', help = 'Please provide a factor to define maxDepth: maxDepth_factor * N = maxDepth. Default = 5')
parser.add_argument('--prefix_extension', '-x', help = 'Option to add to the prefix to distinguish between parameterizations', required = False)

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
		elif ckpt_setting[0] == "qFILTER":
			q = str(ckpt_setting[1])

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

extended_prefix = prefix
writeout_ckpt_filename = args.ckpt_file
if args.prefix_extension is not None:
	extended_prefix = prefix + "-" + args.prefix_extension.replace("_", "-")
	writeout_ckpt_filename = extended_prefix + ".ckpt"
	shutil.copy2(args.ckpt_file, writeout_ckpt_filename)

polymorphic_script = scripts_dir + extended_prefix + "_polymorphicARRAY.sh"
with open(polymorphic_script, 'w') as plm:
	plm.write("#!/bin/bash\n\n")
	plm.write("#SBATCH --cpus-per-task=10\n")
	plm.write("#SBATCH --time=0-20:00:00\n")
	plm.write("#SBATCH --job-name=plm_" + prefix + "\n")
	plm.write("#SBATCH --output=" + jobsout_dir + extended_prefix + "_polymorphic_%A-%a.out\n")
	plm.write("#SBATCH --error=" + jobsout_dir + extended_prefix + "_polymorphic_%A-%a.err\n")
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
		"-out " + gls_dir + extended_prefix + "_${contig}_polymorphic " + \
		"-nThreads 10 " + \
		"-uniqueOnly 1 " + \
		"-remove_bads 1 " + \
		"-trim 0 " + \
		"-C 50 " + \
		"-minMapQ " + args.quality_val + " " + \
		"-minQ " + args.quality_val + " " + \
		"-doCounts 1 " + \
		"-setminDepth " + str(round(float(n_ind) * args.minDepth_factor)) + " " + \
		"-setmaxDepth " + str(round(float(n_ind) * args.maxDepth_factor)) + " " + \
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
	if args.minInd is not None:
		plm.write(" -minInd " + args.minInd)

#Update checkpoint file with the data file names
polymorphic_gls_files = []
polymorphic_maf_files = []
polymorphic_depths_files = []
polymorphic_counts_files = []

for chrom in chrs_list:
	polymorphic_gls_file = gls_dir + extended_prefix + "_" + chrom + "_polymorphic.beagle.gz"
	polymorphic_gls_files.append(polymorphic_gls_file)
	polymorphic_maf_file = gls_dir + extended_prefix + "_" + chrom + "_polymorphic.mafs.gz"
	polymorphic_maf_files.append(polymorphic_maf_file)

if args.global_gls == True:
	global_script = scripts_dir + extended_prefix + "_globalARRAY.sh"
	with open(global_script, 'w') as glb:
		glb.write("#!/bin/bash\n\n")
		glb.write("#SBATCH --cpus-per-task=10\n")
		glb.write("#SBATCH --time=0-20:00:00\n")
		glb.write("#SBATCH --job-name=global_" + prefix + "\n")
		glb.write("#SBATCH --output=" + jobsout_dir + extended_prefix + "_global_%A-%a.out\n")
		glb.write("#SBATCH --error=" + jobsout_dir + extended_prefix + "_global_%A-%a.err\n")
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
			"-out " + gls_dir + extended_prefix + "_${contig}_global " + \
			"-nThreads 10 " + \
			"-uniqueOnly 1 " + \
			"-remove_bads 1 " + \
			"-trim 0 " + \
			"-C 50 " + \
			"-minMapQ " + args.quality_val + " " + \
			"-minQ " + args.quality_val + " " + \
			"-doCounts 1 " + \
			"-setminDepth " + str(round(float(n_ind) * args.minDepth_factor)) + " " + \
			"-setmaxDepth " + str(round(float(n_ind) * args.maxDepth_factor)) + " " + \
			"-GL 1 " + \
			"-doGlf 2 " + \
			"-doMaf 1 " + \
			"-doMajorMinor 1 " + \
			"-doDepth 1 " + \
			"-dumpCounts 3")
		if endedness == "PE":
			glb.write(" -only_proper_pairs 1")
		if args.minInd is not None:
			glb.write(" -minInd " + args.minInd)

	global_gls_files = []
	global_maf_files = []
	global_depths_files = []
	global_counts_files = []

	for chrom in chrs_list:
		global_gls_file = gls_dir + extended_prefix + "_" + chrom + "_global.beagle.gz"
		global_gls_files.append(global_gls_file)
		global_maf_file = gls_dir + extended_prefix + "_" + chrom + "_global.mafs.gz"
		global_maf_files.append(global_maf_file)

with open(writeout_ckpt_filename, 'a') as ckpt:
	ckpt.write("usePREFIX\t" + extended_prefix + "\n")
	ckpt.write("glsDIR\t" + gls_dir + "\n")
	ckpt.write("blackLIST\t" + args.blacklist_inds + "\n")
	ckpt.write("nIND-filtered\t" + str(n_ind) + "\n")
	ckpt.write("bamsLIST-filtered\t" + filtered_bamslist_filename + "\n")
	ckpt.write("depthMIN-factor\t" + str(args.minDepth_factor) + "\n")
	ckpt.write("depthMAX-factor\t" + str(args.maxDepth_factor) + "\n")
	ckpt.write("minQ\t" + args.quality_val + "\n")
	if args.minInd is not None:
		ckpt.write("minIND\t" + args.minInd + "\n")
	ckpt.write("glsFILES-polymorphic\t" + ",".join(polymorphic_gls_files) + "\n")
	ckpt.write("mafsFILES-polymorphic\t" + ",".join(polymorphic_maf_files) + "\n")
	if args.global_gls == True:
		ckpt.write("glsFILES-global\t" + ",".join(global_gls_files) + "\n")
		ckpt.write("mafsFILES-global\t" + ",".join(global_maf_files) + "\n")

print("Step 4 has finished successfully!")
if args.global_gls == True:
	print("You will find two new scripts in ./scripts/:\n" + \
		global_script + " calculates genotype likelihoods across all sites on each chromosome (separately).")
else:
	print("You will find a new script in ./scripts/:")
print(polymorphic_script + " calculates genotype likelihoods across all polymorphic sites on each chromosome (separately).")
print("Both scripts can run simultaneously.")
if args.global_gls == True:
	print("After they have run, you will have genotype likelihoods (gls) and allele frequencies (maf) for " + \
		"all sites in the genome (global) and putatively variable sites (polymorphic).")
else:
	print("After they have run, you will have genotype likelihoods (gls) and allele frequencies (maf) for " + \
		"all putatively variable sites (SNPs).")
if args.prefix_extension is not None:
	print("You have specified a prefix extension. To facilitate downstream branching, a new checkpoint file has been written:")
	print(writeout_ckpt_filename + " contains all the information from " + args.ckpt_file)
print("Additionally, the parameterization used for calculating genotype likelihoods in ANGSD has been recorded in " +\
	writeout_ckpt_filename + ". Including:")
if args.minInd is not None:
	print("the minDepth factor, maxDepth factor, minQ, and minInd")
else:
	print("the minDepth factor, maxDepth factor, and minQ")
