#!/usr/bin/python3

import argparse
from collections import OrderedDict
import os

#Define informative error statements to ensure the config parsing step went according to plan.
def format_path(raw_path, path_type):
	if path_type == "directory":
		if raw_path.startswith("/"):
			if raw_path.endswith("/"):
				formatted_path = raw_path
			else:
				formatted_path = raw_path + "/"
		elif raw_path.startswith("~"):
			if raw_path.endswith("/"):
				formatted_path = os.path.expanduser(raw_path)
			else:
				formatted_path = os.path.expanduser(raw_path) + "/"
	elif path_type == "file":
		if raw_path.startswith("/"):
			formatted_path = raw_path
		elif raw_path.startswith("~"):
			formatted_path = os.path.expanduser(raw_path)
	return formatted_path

def check_wd(wd):
	try:
		assert os.path.isdir(wd), "Path error: " + wd + " is not a valid path. Please provide the path to the directory in which you'd like the results and log files to be saved"
	except AssertionError as path_msg:
		quit(path_msg) #YOURE HERE

def check_input_datafile(user_specified_file):
	try:
		assert os.path.isfile(user_specified_file), "The user-specified file, " + user_specified_file + " does not exist."
	except AssertionError as file_msg:
		quit(file_msg)

def check_fq_readfiles(SEorPE_fastq_files):
	try:
		assert len(SEorPE_fastq_files) < 3, "FASTQs reads error: there are two many files associated with each sample. At most, there should be two files associated with a sample. This may be an issue with parsing the fastq filename: please be sure the raw fastq filenames follow the format <sampleID_R1.fq.gz> or <sampleID_R1.fq> (if not gzipped)."
	except AssertionError as endedness_msg:
		quit(endedness_msg)

def check_fq_filenames(putative_fastq_files):
	for putative_fq in putative_fastq_files:
		try:
			assert os.path.isfile(putative_fq), "FASTQ existential error: the fastq file " + putative_fq + " does not appear to exist."
		except AssertionError as fq_exist_msg:
			quit(fq_exist_msg)
		fq_name_list = putative_fq.split(".")
		try:
			assert fq_name_list[-1] == "fq" or fq_name_list[-2] == "fq" or fq_name_list[-1] == "fastq" or fq_name_list[-2] == "fastq", "FASTQ filename error: the fastq file " + putative_fq + " does not include any of the typical fastq file suffixes (fastq or fq). Please ensure all fastq filenames follow the format <sampleID_R1.fq.gz> or <sampleID_R1.fq> (if not gzipped)."
		except AssertionError as fq_filename_msg:
			quit(fq_filename_msg)

#Read in config file for the run
parser = argparse.ArgumentParser()
parser.add_argument('--config_file', '-c', help = 'Please provide a config file.')
args = parser.parse_args()

#Initialize run config variables with some default values
se_or_pe = None
raw_fastqs_listfile = None
raw_ref_genome = None
raw_working_dir = "~/"
run_prefix = "test"
email = None
adapter_file = None
chrs = "chromosome-1"

#Parse the config file to determine what's needed (if user wants bwa, no need for picard or bowtie2)
with open(args.config_file, 'r') as run_config:
	for raw_config_line in run_config:
		config_line = raw_config_line.rstrip()
		if config_line.startswith('#') == False:
			config_setting = config_line.split('\t')
			if "list of FASTQs" in config_setting[0]:
				raw_fastqs_listfile = config_setting[1]
			elif "file contains adapter sequences for TRIMMOMATIC" in config_setting[0]:
				raw_adapter_file = config_setting[1]
			elif "chromosomes/contigs" in config_setting[0]:
				raw_chrs_file = config_setting[1]
			elif "FASTA file containing the reference genome" in config_setting[0]:
				raw_ref_genome = config_setting[1]
			elif "path to the working directory" in config_setting[0]:
				raw_working_dir = config_setting[1]
			elif "prefix would you like associated with this run" in config_setting[0]:
				run_prefix = config_setting[1].replace("_", "-")
			elif "failed job notifications be sent" in config_setting[0]:
				email = config_setting[1]

#Standardize path formats
my_working_dir = format_path(raw_working_dir, "directory")
fastqs_listfile = format_path(raw_fastqs_listfile, "file")
adapter_file = format_path(raw_adapter_file, "file")
chrs_file = format_path(raw_chrs_file, "file")
ref_genome = format_path(raw_ref_genome, "file")

#Check the parsed info (files, etc) to be sure they are formatted correctly, etc. and throw helpful errors if they aren't.
if check_wd(my_working_dir) is not None:
	print(check_wd(my_working_dir))
if check_input_datafile(fastqs_listfile) is not None:
	print(check_input_datafile(fastqs_listfile))
	print("This file should contain a list of all fastq files that will be assembled.")
if check_input_datafile(adapter_file) is not None:
	print(check_input_datafile(adapter_file))
	print("This file should specify the adapter sequences for TRIMMOMATIC.")
if check_input_datafile(chrs_file) is not None:
	print(check_input_datafile(chrs_file))
	print("This file should specify the chromosomes/scaffolds to be included in the analysis.")
if check_input_datafile(ref_genome) is not None:
	print(check_input_datafile(ref_genome))
	print("This file should specify the reference genome - uncompressed and in FASTA format.")

#Set up a directory structure for results to get printed to
scripts_dir = my_working_dir + "scripts/"
jobsout_dir = my_working_dir + "job_outfiles/"
bwa_dir = my_working_dir + "bwa/"

list_of_dirs = [scripts_dir, jobsout_dir, bwa_dir]
for new_dir in list_of_dirs:
	if os.path.isdir(new_dir) is not True:
		os.mkdir(new_dir)

#Format the chromosomes/scaffolds list for the ckpt file
chrslist = []
with open(chrs_file, 'r') as chromfile:
	for raw_chrom_line in chromfile:
		chrslist.append(raw_chrom_line.rstrip())
chrs = ",".join(chrslist)

#Get a handle on the targeted fastqs
path_to_fastqs_as_list = fastqs_listfile.split("/")
path_to_fastqs = "/".join(path_to_fastqs_as_list[:-1])
fastqs = OrderedDict()
with open(fastqs_listfile, 'r') as listfile:
	for raw_fastq_line in listfile:
		fastqprefix_end = None
		a_fastqfile = raw_fastq_line.rstrip()
		fastqprefix_end = a_fastqfile.split("/")[-1].split("_")
		if fastqprefix_end[0] in fastqs:
			fastqs[fastqprefix_end[0]].append(a_fastqfile)
			se_or_pe = "PE"
		else:
			fastqs[fastqprefix_end[0]] = [a_fastqfile]
			se_or_pe = "SE"

#Index the reference genome
ref_genome_path_list = ref_genome.split("/")
ref_genome_filename = ref_genome_path_list[-1]
ref_genome_filename_as_list = ref_genome_filename.split(".")
refgenome_prefix = ".".join(ref_genome_filename_as_list[:-1])
bwa_script = scripts_dir + refgenome_prefix + "_bwa-indexSLURM.sh"
with open(bwa_script, 'w') as s:
	s.write("#!/bin/bash\n\n")
	s.write("#SBATCH --cpus-per-task=1\n")
	s.write("#SBATCH --job-name=bwa_index_" + refgenome_prefix + "\n")
	s.write("#SBATCH --mail-type=FAIL\n")
	s.write("#SBATCH --mail-user=" + email + "\n")
	s.write("#SBATCH --output=" + jobsout_dir + "bwa-index_" + refgenome_prefix + ".out\n")
	s.write("#SBATCH --error=" + jobsout_dir + "bwa-index_" + refgenome_prefix + ".err\n\n")
	s.write("module unload aligners/bwa/0.7.17\n")
	s.write("module load aligners/bwa/0.7.17\n\n")
	s.write("bwa index -p " + bwa_dir + refgenome_prefix + " " + ref_genome + "\n")

if os.path.isfile(ref_genome_filename + ".fai") == False:
	fai_script = scripts_dir + refgenome_prefix + "_faiSLURM.sh"
	with open(fai_script, 'w') as s:
		s.write("#!/bin/bash\n\n")
		s.write("#SBATCH --cpus-per-task=1\n")
		s.write("#SBATCH --job-name=fai_" + refgenome_prefix + "\n")
		s.write("#SBATCH --mail-type=FAIL\n")
		s.write("#SBATCH --mail-user=" + email + "\n")
		s.write("#SBATCH --output=" + jobsout_dir + "fai_" + refgenome_prefix + ".out\n")
		s.write("#SBATCH --error=" + jobsout_dir + "fai_" + refgenome_prefix + ".err\n")
		s.write("module unload bio/samtools/1.11\n")
		s.write("module load bio/samtools/1.11\n\n")
		s.write("samtools faidx " + ref_genome + "\n")

#Start the checkpoint file and fill it with stuff
with open(run_prefix + ".ckpt", 'w') as ckpt_file:
#capture the directory structure
	ckpt_file.write("workingDIR\t" + my_working_dir + "\n")
	ckpt_file.write("scriptsDIR\t" + scripts_dir + "\n")
	ckpt_file.write("jobsoutDIR\t" + jobsout_dir + "\n")
	ckpt_file.write("adapterFILE\t" + adapter_file + "\n")
	ckpt_file.write("chrsLIST\t" + chrs + "\n")
	ckpt_file.write("bwaDIR\t" + bwa_dir + "\n")
#capture the other odd bits
	ckpt_file.write("prefix\t" + run_prefix + "\n")
	ckpt_file.write("ENDEDNESS\t" + se_or_pe + "\n")
	ckpt_file.write("email\t" + email + "\n")
#capture the files needed to move forward
	ckpt_file.write("refgenomeFASTA\t" + ref_genome + "\n")
	ckpt_file.write("refgenomeIND\t" + refgenome_prefix + "\n")
	ckpt_file.write("nIND\t" + str(len(fastqs.keys())) + "\n")
	for sample_id, fastq_files in fastqs.items():
		if check_fq_readfiles(fastq_files) is not None:
			print(check_fq_readfiles(fastq_files))
		if check_fq_filenames(fastq_files) is not None:
			print(check_fq_filenames(fastq_files))
		else:
			ckpt_file.write("FQ\t" + sample_id + "\t" + " ".join(fastq_files) + "\n")

print("Step 0 has finished successfully! You will find two new scripts in ./scripts/: " + \
	bwa_script + " and " + fai_script + ".")
print("Both scripts can run simultaneously with 'sbatch'.")
print(" However, there is no need to wait for these jobs to finish to move to step 1.")
print("To continue on, call step 1 to generate scripts for fastQC and multiQC.")
print("Remember to pass the newly-made checkpoint (.ckpt) file with the '-p' flag.")
