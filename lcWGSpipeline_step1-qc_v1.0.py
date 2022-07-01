#!/usr/bin/python3

import argparse
import os

#Read in checkpoint file from step 0:
parser = argparse.ArgumentParser()
parser.add_argument('--ckpt_file', '-p', help = 'Please provide the checkpoint file from step 0.')
args = parser.parse_args()

#Initialize run config variables with some default values
working_dir = None
scripts_dir = None
jobsout_dir = None
prefix = None
nind = None
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
		elif ckpt_setting[0] == "jobsoutDIR":
			jobsout_dir = ckpt_setting[1]
		elif ckpt_setting[0] == "prefix":
			prefix = ckpt_setting[1]
		elif ckpt_setting[0] == "nIND":
			nind = ckpt_setting[1]
		elif ckpt_setting[0] == "FQ":
			fq_files = ckpt_setting[2].split(" ")
			fastqs_list.extend(fq_files)

#Set up fastqc/raw directory
fastqc_super_dir = working_dir + "fastqc/"
fastqc_dir = fastqc_super_dir + "raw/"
if os.path.isdir(fastqc_super_dir) is not True:
	os.mkdir(fastqc_super_dir)
if os.path.isdir(fastqc_dir) is not True:
	os.mkdir(fastqc_dir)

#Write an input file for the FASTQC job array
fqc_array_input = scripts_dir + prefix + "-raw_fqcARRAY_input.txt"
iterator = 1
with open(fqc_array_input, 'w') as i:
	for fqfile in fastqs_list:
		i.write(str(iterator) + ":" + fqfile + "\n")
		iterator += 1

#Write the FASTQC job array script
array_script = scripts_dir + prefix + "-raw_fastqcARRAY.sh"
with open(array_script, 'w') as a:
	a.write("#!/bin/bash\n\n")
	a.write("#SBATCH --job-name=fqc_array_" + prefix + "\n")
	a.write("#SBATCH --cpus-per-task=1\n")
	a.write("#SBATCH --output=" + jobsout_dir + prefix + "-raw_fastqc_%A-%a.out" + "\n")
	a.write("#SBATCH --time=3-00:00:00\n")
	a.write("#SBATCH --array=1-" + str(iterator - 1) + "%24\n\n")
	a.write("module unload bio/fastqc/0.11.9\n")
	a.write("module load bio/fastqc/0.11.9\n\n")
	a.write("JOBS_FILE=" + fqc_array_input + "\n")
	a.write("IDS=$(cat ${JOBS_FILE})\n\n")
	a.write("for sample_line in ${IDS}\n")
	a.write("do\n")
	a.write("""\tjob_index=$(echo ${sample_line} | awk -F ":" '{print $1}')\n""")
	a.write("""\tfq=$(echo ${sample_line} | awk -F ":" '{print $2}')\n""")
	a.write("\tif [[ ${SLURM_ARRAY_TASK_ID} == ${job_index} ]]; then\n")
	a.write("\t\tbreak\n")
	a.write("\tfi\n")
	a.write("done\n\n")
	a.write("fastqc ${fq} -o " + fastqc_dir)

mQC_script = scripts_dir + prefix + "-raw_multiqcSLURM.sh"
with open(mQC_script, 'w') as ms:
	ms.write("#!/bin/bash\n\n")
	ms.write("#SBATCH --cpus-per-task=1\n")
	ms.write("#SBATCH --job-name=multiQC\n")
	ms.write("#SBATCH --output=" + jobsout_dir + prefix + "-raw_multiQC.out\n\n")
	ms.write("source /home/ltimm/bin/hydraQC/bin/activate\n")
	ms.write("multiqc " + fastqc_dir)

with open(args.ckpt_file, 'a') as ckpt:
	ckpt.write("fastqc-rawDIR\t" + fastqc_dir + "/n")

print("Step 1 has finished successfully! You will find two new scripts in ./scripts/: " + \
	array_script + " and " + mQC_script + ".")
print(array_script + " must run prior to launching " + mQC_script + ".")
print(" However, there is no need to wait for these jobs to finish to move to step 2.")
print("To continue on, call step 2 to generate a script for TRIMMOMATIC.")
print("Remember to pass the checkpoint (.ckpt) file with the '-p' flag.")
