#!/usr/bin/python3

#does the user need to get the fasta_generate_regions.py script themselves? Or does it come with freebayes?

import argparse
import os
import shutil
from collections import OrderedDict

#Read in config file for the run
parser = argparse.ArgumentParser()
parser.add_argument('--ckpt_file', '-p', help = 'Please provide the checkpoint file created in step 0.')
parser.add_argument('--blacklist_inds', '-b', help = 'Please provide a file listing any individuals that should be removed from downstream analyses.')
parser.add_argument('--prefix_extension', '-x', help = 'Option to add to the prefix to distinguish between parameterizations', required = False)
parser.add_argument('--region_generating_script', '-s', help = 'Provide the location of fasta_generate_regions.py transferred from the freebayes github.')
parser.add_argument('--region_size', '-r', help = 'Option to specify region size', default = "100000")
parser.add_argument('--ncores', '-n', help = 'Option to specify number of cores over which to parallelize', default = "20")
args = parser.parse_args()

#Initialize run config variables with some default values
working_dir = None
scripts_dir = None
jobsout_dir = None
ref_genome = None
prefix = None
email = None
n_ind = None
blacklist_bams = []
full_bams_list = None

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
		elif ckpt_setting[0] == "nIND":
			n_ind = int(ckpt_setting[1])
		elif ckpt_setting[0] == "bamsLIST-all":
			full_bams_list = ckpt_setting[1]

vcf_dir = working_dir + "vcf/"
if os.path.isdir(vcf_dir) is not True:
	os.mkdir(vcf_dir)

#remove blacklisted samples
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

#copy to new ckpt, if applicable
extended_prefix = prefix
writeout_ckpt_filename = args.ckpt_file
if args.prefix_extension is not None:
	extended_prefix = prefix + "-" + args.prefix_extension.replace("_", "-")
	writeout_ckpt_filename = extended_prefix + ".ckpt"
	shutil.copy2(args.ckpt_file, writeout_ckpt_filename)

#FREEBAYES-PARALLEL
freebayes_parallel_script = scripts_dir + extended_prefix + "_freebayes.sh"
with open(freebayes_parallel_script, 'w') as fbp:
	fbp.write("#!/bin/bash\n\n")
	fbp.write("#SBATCH --cpus-per-task=" + args.ncores + "\n")
	fbp.write("#SBATCH --time=3-00:00:00\n")
	fbp.write("#SBATCH --job-name=fbp_" + extended_prefix + "\n")
	fbp.write("#SBATCH --output=" + jobsout_dir + extended_prefix + "_polymorphic_%A-%a.out\n")
	fbp.write("#SBATCH --error=" + jobsout_dir + extended_prefix + "_polymorphic_%A-%a.err\n")
	fbp.write("#SBATCH --mail-type=FAIL\n")
	fbp.write("#SBATCH --mail-user=" + email + "\n\n")

	fbp.write("mamba activate freebayes-1.3.7\n\n")

	fbp.write("freebayes-parallel <(" + \
		args.region_generating_script + " " + \
		ref_genome + ".fai " + \
		args.region_size + ") " + \
		args.ncores + " " +
		"-f " + ref_genome + " " + \
		"-L " + filtered_bamslist_filename + " > " + \
		vcf_dir + extended_prefix + ".vcf")

#freebayes-parallel <(fasta_generate_regions.py ref.fa.fai 100000) NCORES \
#    -f ref.fa -L bamslist.txt > prefix.vcf
#Update checkpoint file

with open(writeout_ckpt_filename, 'a') as ckpt:
	ckpt.write("vcfDIR\t" + vcf_dir + "\n")
	ckpt.write("blackLIST\t" + args.blacklist_inds + "\n")
	ckpt.write("nIND-filtered\t" + str(n_ind) + "\n")
	ckpt.write("bamsLIST-filtered\t" + filtered_bamslist_filename + "\n")
	ckpt.write("vcfOUT\t" + extended_prefix + ".vcf\n")
	ckpt.write("usePREFIX\t" + extended_prefix + "\n")

print("The freebayes version of Step 4 has finished successfully!")
print(freebayes_parallel_script + " will detect genetic variants and report them in " + extended_prefix + ".vcf\n")