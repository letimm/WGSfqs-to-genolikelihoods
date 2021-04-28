#!/usr/bin/python3

import argparse
import subprocess

#Read in checkpoint file from step 0:
parser = argparse.ArgumentParser()
parser.add_argument('--ckpt_file', '-p', help = 'Please provide the checkpoint file from step 0.')
args = parser.parse_args()

#Initialize run config variables with some default values
working_dir = None
scripts_dir = None
fastqc_dir = None
prefix = None
fastqs_list = []

#Parse the checkpoint file from step 0
with open(args.ckpt_file, 'r') as last_step_ckpt:
	for raw_ckpt_line in last_step_ckpt:
		ckpt_line = raw_ckpt_line.rstrip()
		ckpt_setting = ckpt_line.split('\t')
		if ckpt_setting[0] == "workingDIR":
			working_dir = ckpt_setting[1]
		elif ckpt_setting[0] == "scriptsDIR":
			scripts_dir = ckpt_setting[1]
		elif ckpt_setting[0] == "fastqcDIR":
			fastqc_dir = ckpt_setting[1]
		elif ckpt_setting[0] == "prefix":
			prefix = ckpt_setting[1]
		elif ckpt_setting[0] == "FQ":
			fq_files = ckpt_setting[2].split(" ")
			fastqs_list.extend(fq_files)

#FASTQC all fastqs
for fq in fastqs_list:
	fq_name_as_list = fq.split(".")
	fq_label = fq_name_as_list[0]
	script = scripts_dir + prefix + "_" + fq_label + "_fastqcSLURM.sh"
	with open(script, 'w') as s:
		s.write("#!/bin/bash\n\n")
		s.write("#SBATCH --nodes=1\n")
		s.write("#SBATCH --ntasks=1\n")
		s.write("#SBATCH --job-name=fqc_" + fq_label + "\n")
		s.write("#SBATCH --output=" + fastqc_dir + prefix + "_" + fq_label + "_fastqc.out\n\n")
		s.write("module unload bio/fastqc\n")
		s.write("module load bio/fastqc\n\n")
		s.write("fastqc " + working_dir + fq + " -o " + fastqc_dir + "\n")
	subprocess.call(["sbatch", script])