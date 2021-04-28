#!/usr/bin/python3
import argparse
from collections import OrderedDict
import math
import subprocess

#Read in config file for the run
parser = argparse.ArgumentParser()
parser.add_argument('--ckpt_file', '-p', help = 'Please provide the checkpoint file created in step 0.')
args = parser.parse_args()

#Initialize run config variables with some default values
working_dir = None
scripts_dir = None
bwa_dir = None
samtools_dir = None
max_nodes = None
prefix = None
refgenome_prefix = None
fastqs = OrderedDict()
samfiles = []
endedness = None

#Parse the config file to determine what's needed (if user wants bwa, no need for picard or bowtie2)
with open(args.ckpt_file, 'r') as last_step_ckpt:
	for raw_ckpt_line in last_step_ckpt:
		ckpt_line = raw_ckpt_line.rstrip()
		ckpt_setting = ckpt_line.split('\t')
		if ckpt_setting[0] == "workingDIR":
			working_dir = ckpt_setting[1]
		elif ckpt_setting[0] == "scriptsDIR":
			scripts_dir = ckpt_setting[1]
		elif ckpt_setting[0] == "bwaDIR":
			bwa_dir = ckpt_setting[1]
		elif ckpt_setting[0] == "samtoolsDIR":
			samtools_dir = ckpt_setting[1]
		elif ckpt_setting[0] == "maxNODES":
			max_nodes = float(ckpt_setting[1])
		elif ckpt_setting[0] == "prefix":
			prefix = ckpt_setting[1]
		elif ckpt_setting[0] == "refgenomeIND":
			refgenome_prefix = ckpt_setting[1]
		elif ckpt_setting[0] == "FQ":
			fastqs[ckpt_setting[1]] = ckpt_setting[2].split(" ")

#Determine how many nodes can be used to run each job
n_nodes = math.floor(max_nodes / float(len(fastqs)))
if max_nodes == 1:
	n_nodes = 1
elif n_nodes <= 2:
	n_nodes = 2
else:
	n_nodes = n_nodes

#Align each sample's fastq(s) to the reference genome
for sample_id, fq_files in fastqs.items():
	samfile = prefix + "_" + sample_id + ".sam"
	script = scripts_dir + prefix + "_" + sample_id + "_bwa-memSLURM.sh"
	samfiles.append(samfile)
	with open(script, 'w') as s:
		s.write("#!/bin/bash\n\n")
		s.write("#SBATCH --nodes=" + str(n_nodes) + "\n")
		s.write("#SBATCH --ntasks=" + str(n_nodes) + "\n")
		s.write("#SBATCH --job-name=bwa-mem_" + sample_id + "\n")
		s.write("#SBATCH --output=" + bwa_dir + prefix + "_" + sample_id + "_bwa-mem.err\n\n")
		s.write("module unload aligners/bwa\n")
		s.write("module load aligners/bwa\n\n")

		if len(fq_files) == 1:
			endedness = "SE"
			s.write("bwa mem -M -t 1 " + bwa_dir + refgenome_prefix + " " + working_dir + fq_files[0] + " 2> " + \
				bwa_dir + prefix + "_" + sample_id + "_bwa-mem.out > " + samtools_dir + samfile + "\n")
		elif len(fq_files) == 2:
			endedness = "PE"
			s.write("bwa mem -M -t 1 " + bwa_dir + refgenome_prefix + " " + working_dir + fq_files[0] + working_dir + fq_files[1] + " 2> " + \
				bwa_dir + prefix + "_" + sample_id + "_bwa-mem.out > " + samtools_dir + samfile + "\n")
		else:
			print("Whoops. There are either >2 fq files associated with a sample or there are no fq files associated with a sample. Look into that, please.")
	subprocess.call(["sbatch", script])

#Update checkpoint file with .sam filenames resulting from alignment with bwa
with open(args.ckpt_file, 'a') as ckpt:
	ckpt.write("distNODES\t" + str(n_nodes) + "\n")
	ckpt.write("ENDEDNESS\t" + endedness + "\n")
	for sam in samfiles:
		ckpt.write("SAM\t" + sam + "\n")
