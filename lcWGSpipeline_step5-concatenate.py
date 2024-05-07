#!/usr/bin/python3
import argparse

#Read in config file for the run
parser = argparse.ArgumentParser()
parser.add_argument('--ckpt_file', '-p', help = 'Please provide the checkpoint file created in step 0.')
args = parser.parse_args()

scripts_dir = None
jobsout_dir = None
prefix = None
email = None
gls_dir = None
polymorphic_beagle_files = []
polymorphic_maf_files = []
wgp_beagle = None
wgp_maf = None

#Parse the config file to determine what's needed
with open(args.ckpt_file, 'r') as last_step_ckpt:
	for raw_ckpt_line in last_step_ckpt:
		ckpt_line = raw_ckpt_line.rstrip()
		ckpt_setting = ckpt_line.split('\t')
		if ckpt_setting[0] == "scriptsDIR":
			scripts_dir = ckpt_setting[1]
		elif ckpt_setting[0] == "jobsoutDIR":
			jobsout_dir = ckpt_setting[1]
		elif ckpt_setting[0] == "glsDIR":
			gls_dir = ckpt_setting[1]
		elif ckpt_setting[0] == "usePREFIX":
			prefix = ckpt_setting[1]
		elif ckpt_setting[0] == "email":
			email = ckpt_setting[1]
		elif ckpt_setting[0] == "glsFILES-polymorphic":
			polymorphic_beagle_files = ckpt_setting[1].split(",")
		elif ckpt_setting[0] == "mafsFILES-polymorphic":
			polymorphic_maf_files = ckpt_setting[1].split(",")

#Concatenate beagle files to make the polymorphic gls file for the whole genome
all_chrs_polymorphic_beagle = gls_dir + prefix + "_wholegenome_polymorphic.beagle"
beagle_concat_script = scripts_dir + prefix + "_concatenate_beagles.sh"
with open(beagle_concat_script, 'w') as bcs:
	bcs.write("#!/bin/bash\n\n")
	bcs.write("#SBATCH --cpus-per-task=10\n")
	bcs.write("#SBATCH --job-name=" + prefix + "_concat-beagles\n")
	bcs.write("#SBATCH --mail-type=FAIL\n")
	bcs.write("#SBATCH --mail-user=" + email + "\n")
	bcs.write("#SBATCH --output=" + jobsout_dir + prefix + "_concatenate-beagles_%A.out\n")
	bcs.write("#SBATCH --error=" + jobsout_dir + prefix + "_concatenate-beagles_%A.err\n\n")
	bcs.write("zcat " + polymorphic_beagle_files[0] + " " + \
		"| head -n 1 > " + \
		all_chrs_polymorphic_beagle + "; " + \
		"for i in " + " ".join(polymorphic_beagle_files) + "; " + \
		"do zcat $i " + \
		"| tail -n +2 -q >> " + \
		all_chrs_polymorphic_beagle + "; " + \
		"done\n")
	bcs.write("gzip " + all_chrs_polymorphic_beagle)

#Concatenate maf files to make the polymorphic maf file for the whole genome
all_chrs_polymorphic_maf = gls_dir + prefix + "_wholegenome_polymorphic.mafs"
all_chrs_polymorphic_sites = gls_dir + prefix + "_wholegenome_polymorphic.sites"
mafs_concat_script = scripts_dir + prefix + "_concatenate_mafs.sh"
with open(mafs_concat_script, 'w') as mcs:
	mcs.write("#!/bin/bash\n\n")
	mcs.write("#SBATCH --cpus-per-task=10\n")
	mcs.write("#SBATCH --job-name=" + prefix + "_concat-mafs\n")
	mcs.write("#SBATCH --mail-type=FAIL\n")
	mcs.write("#SBATCH --mail-user=" + email + "\n")
	mcs.write("#SBATCH --output=" + jobsout_dir + prefix + "_concatenate-mafs_%A.out\n")
	mcs.write("#SBATCH --error=" + jobsout_dir + prefix + "_concatenate-mafs_%A.err\n\n")
	mcs.write("module unload bio/angsd/0.933\n")
	mcs.write("module load bio/angsd/0.933\n\n")
	mcs.write("for i in " + " ".join(polymorphic_maf_files) + "\n" + \
		#this is a little dangerous - if a mafs file already exists, it will just append to the existing and you'll end up with 2x, 3x, etc sites
		"do zcat $i | tail -n +2 -q >> " + \
		all_chrs_polymorphic_maf + "; " + \
		"done\n")
	mcs.write("cut -f 1,2,3,4 " + all_chrs_polymorphic_maf + " > " + all_chrs_polymorphic_sites + "\n")
	mcs.write("gzip " + all_chrs_polymorphic_maf + "\n\n")
	mcs.write("angsd sites index " + all_chrs_polymorphic_sites + "\n")

with open(args.ckpt_file, 'a') as ckpt:
	ckpt.write("wholegenomeBEAGLE-polymorphic\t" + all_chrs_polymorphic_beagle + ".gz\n")
	ckpt.write("wholegenomeMAF-polymorphic\t" + all_chrs_polymorphic_maf + ".gz\n")
	ckpt.write("wholegenomeSITES-polymorphic\t" + all_chrs_polymorphic_sites + "\n")

print("Two new scripts have been generated: \n")
print(beagle_concat_script + " concatenates all the polymorphic beagle files to produce a single file containing genotype likelihoods for all polymorphic sites across the genome.\n")
print(mafs_concat_script + " concatenates all the polymorphic maf files to produce a single file containing minor allele frequencies for all polymorphic sites across the genome.\n")
print("Congratulations! Data assembly is complete!")
