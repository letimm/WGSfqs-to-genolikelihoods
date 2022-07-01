#!/usr/bin/python3

import argparse
import os
from collections import OrderedDict
import math

#Read in config file for the run
parser = argparse.ArgumentParser()
parser.add_argument('--ckpt_file', '-p', help = 'Please provide the checkpoint file created in step 0.')
args = parser.parse_args()

#Initialize run config variables with some default values
working_dir = None
scripts_dir = None
jobsout_dir = None
endedness = None
adapter_file = None
prefix = None
fastqs = OrderedDict()
trimmed_files = OrderedDict()

#Parse the config file to determine what's needed (if user wants bwa, no need for picard or bowtie2)
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
		elif ckpt_setting[0] == "ENDEDNESS":
			endedness = ckpt_setting[1]
		elif ckpt_setting[0] == "adapterFILE":
			adapter_file = ckpt_setting[1]
		elif ckpt_setting[0] == "prefix":
			prefix = ckpt_setting[1]
		elif ckpt_setting[0] == "FQ":
			fastqs[ckpt_setting[1]] = ckpt_setting[2].split(" ")

#Set up fastqc/trimmed directory
fastqc_super_dir = working_dir + "fastqc/"
fastqc_dir = fastqc_super_dir + "trimmed/"
if os.path.isdir(fastqc_super_dir) is not True:
	os.mkdir(fastqc_super_dir)
if os.path.isdir(fastqc_dir) is not True:
	os.mkdir(fastqc_dir)

##Set up trimmed directory directory
trim_dir = working_dir + "trimmed/"
if os.path.isdir(trim_dir) is not True:
	os.mkdir(trim_dir)

#Write an input file for the alignment job array
trim_array_input = scripts_dir + prefix + "_trimARRAY_input.txt"
with open(trim_array_input, 'w') as i:
	if endedness == "SE":
		iterator = 1
		for se_fastq in fastqs.values():
			i.write(str(iterator) + ":" + se_fastq + "\n")
			iterator += 1
	elif endedness == "PE":
		iterator = 1
		for pe_fastqs in fastqs.values():
			i.write(str(iterator) + ":" + ":".join(pe_fastqs) + "\n")
			iterator += 1

#Write the trimmomatic job array script
trim_array_script = scripts_dir + prefix + "_trimARRAY.sh"
file_suffix = None
with open(trim_array_script, 'w') as t:
	t.write("#!/bin/bash\n\n")
	t.write("#SBATCH --job-name=trim\n")
	t.write("#SBATCH --cpus-per-task=4\n")
	t.write("#SBATCH --output=" + jobsout_dir + prefix + "_trimming_%A-%a.out" + "\n")
	t.write("#SBATCH --time=7-00:00:00\n")
	t.write("#SBATCH --array=1-" + str(len(fastqs.keys())) + "%48\n\n")
	t.write("module unload bio/trimmomatic/0.39\n")
	t.write("module load bio/trimmomatic/0.39\n\n")
	t.write("JOBS_FILE=" + trim_array_input + "\n")
	t.write("IDS=$(cat ${JOBS_FILE})\n\n")
	t.write("for sample_line in ${IDS}\n")
	t.write("do\n")
	t.write("""\tjob_index=$(echo ${sample_line} | awk -F ":" '{print $1}')\n""")
	t.write("""\tfq_r1=$(echo ${sample_line} | awk -F ":" '{print $2}')\n""")
	if endedness == "PE":
		t.write("""\tfq_r2=$(echo ${sample_line} | awk -F ":" '{print $3}')\n""")
	t.write("\tif [[ ${SLURM_ARRAY_TASK_ID} == ${job_index} ]]; then\n")
	t.write("\t\tbreak\n")
	t.write("\tfi\n")
	t.write("done\n\n")

	t.write("sample_id=$(echo $fq_r1 | sed 's!^.*/!!')\n")
	t.write("sample_id=${sample_id%%_*}\n\n")
	if endedness == "SE":
		t.write("java -jar ${TRIMMOMATIC} " + \
			"SE " + \
			"-threads 4 " + \
			"-phred33 " + \
			"${fq_r1} " + \
			trim_dir + "${sample_id}_trimmed.fq.gz " + \
			"ILLUMINACLIP:" + adapter_file + ":2:30:10 MINLEN:40\n")
	elif endedness == "PE":
		t.write("java -jar ${TRIMMOMATIC} " + \
			"PE " + \
			"-threads 4 " + \
			"-phred33 " + \
			"${fq_r1} " + \
			"${fq_r2} " + \
			trim_dir + "${sample_id}_trimmed_R1_paired.fq.gz " + \
			trim_dir + "${sample_id}_trimmed_R1_unpaired.fq.gz " + \
			trim_dir + "${sample_id}_trimmed_R2_paired.fq.gz " + \
			trim_dir + "${sample_id}_trimmed_R2_unpaired.fq.gz " + \
			"ILLUMINACLIP:" + adapter_file + ":2:30:10:1:true MINLEN:40\n")

#Rerun FASTQC and MULTIQC now that TRIMMOMATIC has run
fastqc_array_input = scripts_dir + prefix + "-trim_fqcARRAY_input.txt"
iterator = 1
with open(fastqc_array_input, 'w') as i:
	for fastq in fastqs.keys():
		i.write(str(iterator) + ":" + trim_dir + fastq + "_trimmed_R1_paired.fq.gz\n")
		iterator += 1
		if endedness == "PE":
			i.write(str(iterator) + ":" + trim_dir + fastq + "_trimmed_R2_paired.fq.gz\n")
			iterator += 1

#Write the depth calculations job array script
fq_array_script = scripts_dir + prefix + "-trim_fastqcARRAY.sh"
with open(fq_array_script, 'w') as f:
	f.write("#!/bin/bash\n\n")
	f.write("#SBATCH --job-name=fqc_array_" + prefix + "\n")
	f.write("#SBATCH --cpus-per-task=1\n")
	f.write("#SBATCH --output=" + jobsout_dir + prefix + "-trim_fastqc_%A-%a.out" + "\n")
	f.write("#SBATCH --time=3-00:00:00\n")
	f.write("#SBATCH --array=1-" + str(iterator - 1) + "%24\n\n")
	f.write("module unload bio/fastqc/0.11.9\n")
	f.write("module load bio/fastqc/0.11.9\n\n")
	f.write("JOBS_FILE=" + fastqc_array_input + "\n")
	f.write("IDS=$(cat ${JOBS_FILE})\n\n")
	f.write("for sample_line in ${IDS}\n")
	f.write("do\n")
	f.write("""\tjob_index=$(echo ${sample_line} | awk -F ":" '{print $1}')\n""")
	f.write("""\tfq=$(echo ${sample_line} | awk -F ":" '{print $2}')\n""")
	f.write("\tif [[ ${SLURM_ARRAY_TASK_ID} == ${job_index} ]]; then\n")
	f.write("\t\tbreak\n")
	f.write("\tfi\n")
	f.write("done\n\n")
	f.write("fastqc ${fq} -o " + fastqc_dir)

mQC_script = scripts_dir + prefix + "-trim_multiqcSLURM.sh"
with open(mQC_script, 'w') as ms:
	ms.write("#!/bin/bash\n\n")
	ms.write("#SBATCH --cpus-per-task=1\n")
	ms.write("#SBATCH --job-name=multiQC\n")
	ms.write("#SBATCH --output=" + jobsout_dir + prefix + "-trim_multiQC.out\n\n")
	ms.write("source /home/ltimm/bin/hydraQC/bin/activate\n")
	ms.write("multiqc " + fastqc_dir)

#Update checkpoint file with trimmed filenames resulting from alignment
with open(args.ckpt_file, 'a') as ckpt:
	ckpt.write("fastqc-trimDIR\t" + fastqc_dir + "\n")
	ckpt.write("trimmedDIR\t" + trim_dir + "\n")
	if endedness == "SE":
		for fq_id in fastqs.keys():
			ckpt.write("trimmedFQ\t" + fq_id + \
				"\t" + trim_dir + fq_id + "_trimmed.fq.gz\n")
	elif endedness == "PE":
		for fq_id in fastqs.keys():
			ckpt.write("trimmedFQ\t" + fq_id + \
				"\t" + trim_dir + fq_id + "_trimmed_R1_paired.fq.gz " + \
				trim_dir + fq_id + "_trimmed_R2_paired.fq.gz\n")

print("Step 2 has finished successfully! You will find three new scripts in ./scripts/: " + \
	trim_array_script + ", " + fq_array_script + ", and " + mQC_script + ".")
print("Each script must run in the above order and each job must finish before submitting the next.")
print("While there is no need to wait before running step 3 to generate the alignment script, " + \
	"it is probably wise to wait until the multiQC script has completed and the results have been viewed " + \
	"in a web browser prior to submitting the script written by step 3.")
print("Remember to pass the checkpoint (.ckpt) file with the '-p' flag.")
