#!/usr/bin/python3
#let's start this as a python thing; if bash has more traction, we'll do it in bash! (or GOlang)

import argparse
from collections import OrderedDict
import os
import subprocess

#Define informative error statements to ensure the config parsing step went according to plan.
def check_cluster_resources(nnodes):
	try:
		assert nnodes <= 128
	except ValueError:
		return "Nodes error: you have specified >128 nodes, which is too many nodes per core. Please specify 128 or fewer nodes."

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
		assert os.path.isdir(wd)
	except ValueError:
		return "Path error: " + wd + " is not a valid path. Please provide the path to the directory in which you'd like the results and logfiles to be saved"

def check_input_datafiles(fastqslist, refgenfasta):
	try:
		assert os.path.isfile(fastqslist)
	except ValueError:
		return "FASTQs list file error: the file containing a list of FASTQs does not exist. Please specify a file with the FASTQs to be analyzed."
	try:
		assert os.path.isfile(refgenfasta)
	except ValueError:
		return "Reference genome file error: the reference genome FASTA file does not exist. Please specify a FASTA file with the reference genome."

def check_fq_readfiles(SEorPE_fastq_files):
	try:
		assert len(SEorPE_fastq_files) < 3
	except ValueError:
		return "FASTQs reads error: there are two many files associated with each sample. At most, there should be two files associated with a sample. This may be an issue with parsing the fastq filename: please be sure the raw fastq filenames follow the format <sampleID_R1.fq.gz> or <sampleID_R1.fq> (if not gzipped)."

def check_fq_filenames(putative_fastq_files):
	for putative_fq in putative_fastq_files:
		try:
			assert os.path.isfile(putative_fq)
		except ValueError:
			return "FASTQ existential error: the fastq file " + putative_fq + " does not appear to exist."
		fq_name_list = putative_fq.split(".")
		try:
			assert fq_name_list[-1] == "fq" or fq_name_list[-2] == "fq" or fq_name_list[-1] == "fastq" or fq_name_list[-2] == "fastq"
		except ValueError:
			return "FASTQ filename error: the fastq file " + putative_fq + " does not include any of the typical fastq file suffixes (fastq or fq). Please ensure all fastq filenames follow the format <sampleID_R1.fq.gz> or <sampleID_R1.fq> (if not gzipped)."

#Read in config file for the run
parser = argparse.ArgumentParser()
parser.add_argument('--config_file', '-c', help = 'Please provide a config file.')
args = parser.parse_args()

#Initialize run config variables with some default values
max_nodes = "1"
se_or_pe = "se"
raw_fastqs_listfile = ""
raw_ref_genome = ""
raw_working_dir = "~/"
run_prefix = "test"

#Parse the config file to determine what's needed (if user wants bwa, no need for picard or bowtie2)
with open(args.config_file, 'r') as run_config:
	for raw_config_line in run_config:
		config_line = raw_config_line.rstrip()
		if config_line.startswith('#') == False:
			config_setting = config_line.split('\t')
			if "nodes can you use (max)" in config_setting[0]:
				max_nodes = int(config_setting[1])
			elif "list of FASTQs" in config_setting[0]:
				raw_fastqs_listfile = config_setting[1]
			elif "FASTA file containing the reference genome" in config_setting[0]:
				raw_ref_genome = config_setting[1]
			elif "path to the working directory" in config_setting[0]:
				raw_working_dir = config_setting[1]
			elif "prefix would you like associated with this run" in config_setting[0]:
				run_prefix = config_setting[1]
			#Will probably have to add more or rearrange; but this is good for now.

#Standardize path formats
my_working_dir = format_path(raw_working_dir, "directory")
fastqs_listfile = format_path(raw_fastqs_listfile, "file")
ref_genome = format_path(raw_ref_genome, "file")

#Check the parsed info (files, etc) to be sure they are formatted correctly, etc. and throw helpful errors if they aren't.
if check_cluster_resources(max_nodes) is not None:
	print(check_cluster_resources(max_nodes))
if check_wd(my_working_dir) is not None:
	print(check_wd(my_working_dir))
if check_input_datafiles(fastqs_listfile, ref_genome) is not None:
	print(check_input_datafiles(fastqs_listfile, ref_genome))

#Set up a directory structure for results to get printed to
scripts_dir = my_working_dir + "scripts/"
fastqc_dir = my_working_dir + "fastqc/"
bwa_dir = my_working_dir + "bwa/"
samtools_dir = my_working_dir + "samtools/"
bamtools_dir = my_working_dir + "bamtools/"
angsd_dir = my_working_dir + "angsd/"

if os.path.isdir(scripts_dir) is not True:
	os.mkdir(scripts_dir)
if os.path.isdir(fastqc_dir) is not True:
	os.mkdir(fastqc_dir)
if os.path.isdir(bwa_dir) is not True:
	os.mkdir(bwa_dir)
if os.path.isdir(samtools_dir) is not True:
	os.mkdir(samtools_dir)
if os.path.isdir(bamtools_dir) is not True:
	os.mkdir(bamtools_dir)
if os.path.isdir(angsd_dir) is not True:
	os.mkdir(angsd_dir)

#Get a handle on the targeted fastqs
path_to_fastqs_as_list = fastqs_listfile.split("/")
path_to_fastqs = "/".join(path_to_fastqs_as_list[:-1])
fastqs = OrderedDict()
with open(fastqs_listfile, 'r') as listfile:
	for raw_fastq_line in listfile:
		a_fastqfile = raw_fastq_line.rstrip()
		fastqprefix_read = a_fastqfile.split("_")
		if fastqprefix_read[0] in fastqs:
			fastqs[fastqprefix_read[0]].append(a_fastqfile)
		else:
			fastqs[fastqprefix_read[0]] = [a_fastqfile]

#Index the reference genome
ref_genome_path_list = ref_genome.split("/")
ref_genome_filename = ref_genome_path_list[-1]
ref_genome_filename_as_list = ref_genome_filename.split(".")
refgenome_prefix = ".".join(ref_genome_filename_as_list[:-1])
bwa_script = scripts_dir + refgenome_prefix + "_bwa-indexSLURM.sh"
with open(bwa_script, 'w') as s:
	s.write("#!/bin/bash\n\n")
	s.write("#SBATCH --nodes=1\n")
	s.write("#SBATCH --ntasks=1\n")
	s.write("#SBATCH --job-name=bwa_index_" + refgenome_prefix + "\n")
	s.write("#SBATCH --output=" + bwa_dir + "bwa-index_" + refgenome_prefix + ".out\n\n")
	s.write("module unload aligners/bwa\n")
	s.write("module load aligners/bwa\n\n")
	s.write("bwa index -p " + bwa_dir + refgenome_prefix + " " + ref_genome + "\n")
subprocess.call(["sbatch", bwa_script])

if os.path.isfile(ref_genome_filename + ".fai") == False:
	fai_script = scripts_dir + refgenome_prefix + "_faiSLURM.sh"
	with open(fai_script, 'w') as s:
		s.write("#!/bin/bash\n\n")
		s.write("#SBATCH --nodes=1\n")
		s.write("#SBATCH --ntasks=1\n")
		s.write("#SBATCH --job-name=fai_" + refgenome_prefix + "\n")
		s.write("#SBATCH --output=" + samtools_dir + "fai_" + refgenome_prefix + ".out\n\n")
		s.write("module unload bio/samtools\n")
		s.write("module load bio/samtools\n\n")
		s.write("samtools faidx " + ref_genome + "\n")
	subprocess.call(["sbatch", fai_script])

#Start the checkpoint file and fill it with stuff from step0
with open(run_prefix + ".ckpt", 'w') as ckpt_file:
#capture the directory structure
	ckpt_file.write("workingDIR\t" + my_working_dir + "\n")
	ckpt_file.write("scriptsDIR\t" + scripts_dir + "\n")
	ckpt_file.write("fastqcDIR\t" + fastqc_dir + "\n")
	ckpt_file.write("bwaDIR\t" + bwa_dir + "\n")
	ckpt_file.write("samtoolsDIR\t" + samtools_dir + "\n")
	ckpt_file.write("bamtoolsDIR\t" + bamtools_dir + "\n")
	ckpt_file.write("angsdDIR\t" + angsd_dir + "\n")
#capture the other odd bits
	ckpt_file.write("maxNODES\t" + str(max_nodes) + "\n")
	ckpt_file.write("prefix\t" + run_prefix + "\n")
#capture the files needed to move forward
	ckpt_file.write("refgenomeFASTA\t" + ref_genome + "\n")
	ckpt_file.write("refgenomeIND\t" + refgenome_prefix + "\n")
	for sample_id, fastq_files in fastqs.items():
		if check_fq_readfiles(fastq_files) is not None:
			print(check_fq_readfiles(fastq_files))
		if check_fq_filenames(fastq_files) is not None:
			print(check_fq_filenames(fastq_files))
		ckpt_file.write("FQ\t" + sample_id + "\t" + " ".join(fastq_files) + "\n")
