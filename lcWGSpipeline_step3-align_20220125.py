#!/usr/bin/python3

import argparse
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
bwa_dir = None
samtools_dir = None
bamtools_dir = None
prefix = None
refgenome_prefix = None
trimmed_fastqs = OrderedDict()
endedness = None
depth_files = []

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
		elif ckpt_setting[0] == "bwaDIR":
			bwa_dir = ckpt_setting[1]
		elif ckpt_setting[0] == "samtoolsDIR":
			samtools_dir = ckpt_setting[1]
		elif ckpt_setting[0] == "bamtoolsDIR":
			bamtools_dir = ckpt_setting[1]
		elif ckpt_setting[0] == "prefix":
			prefix = ckpt_setting[1]
		elif ckpt_setting[0] == "ENDEDNESS":
			endedness = ckpt_setting[1]
		elif ckpt_setting[0] == "refgenomeIND":
			refgenome_prefix = ckpt_setting[1]
		elif ckpt_setting[0] == "trimmedFQ":
			trimmed_fastqs[ckpt_setting[1]] = ckpt_setting[2].split(" ")

#Write an input file for the alignment job array
align_array_input = scripts_dir + prefix + "_alignARRAY_input.txt"
with open(align_array_input, 'w') as i:
	if endedness == "SE":
		iterator = 1
		for se_fastq in trimmed_fastqs.values():
			i.write(str(iterator) + ":" + se_fastq + "\n")
			iterator += 1
	elif endedness == "PE":
		iterator = 1
		for pe_fastqs in trimmed_fastqs.values():
			i.write(str(iterator) + ":" + ":".join(pe_fastqs) + "\n")
			iterator += 1

#Write the alignment job array script
align_array_script = scripts_dir + prefix + "_alignARRAY.sh"
file_suffix = None
with open(align_array_script, 'w') as a:
	a.write("#!/bin/bash\n\n")
	a.write("#SBATCH --job-name=align\n")
	a.write("#SBATCH --cpus-per-task=10\n")
	a.write("#SBATCH --output=" + jobsout_dir + prefix + "_alignment_%A-%a.out" + "\n")
	a.write("#SBATCH --time=7-00:00:00\n")
	a.write("#SBATCH --array=1-" + str(len(trimmed_fastqs.keys())) + "%12\n\n")
	a.write("module unload aligners/bwa/0.7.17 bio/samtools/1.11 bio/bamtools/2.5.1 bio/picard/2.23.9 bio/bamutil/1.0.5\n")
	a.write("module load aligners/bwa/0.7.17 bio/samtools/1.11 bio/bamtools/2.5.1 bio/picard/2.23.9 bio/bamutil/1.0.5\n\n")
	a.write("JOBS_FILE=" + align_array_input + "\n")
	a.write("IDS=$(cat ${JOBS_FILE})\n\n")
	a.write("for sample_line in ${IDS}\n")
	a.write("do\n")
	a.write("""\tjob_index=$(echo ${sample_line} | awk -F ":" '{print $1}')\n""")
	a.write("""\tfq_r1=$(echo ${sample_line} | awk -F ":" '{print $2}')\n""")
	if endedness == "PE":
		a.write("""\tfq_r2=$(echo ${sample_line} | awk -F ":" '{print $3}')\n""")
	a.write("\tif [[ ${SLURM_ARRAY_TASK_ID} == ${job_index} ]]; then\n")
	a.write("\t\tbreak\n")
	a.write("\tfi\n")
	a.write("done\n\n")

	a.write("sample_id=$(echo $fq_r1 | sed 's!^.*/!!')\n")
	a.write("sample_id=${sample_id%%_*}\n\n")

	if endedness == "SE":
		a.write("bwa mem -M -t 10 " + bwa_dir + refgenome_prefix + " ${fq_r1} 2> " + \
				bwa_dir + prefix + "_${sample_id}_bwa-mem.out > " + \
				samtools_dir + prefix + "_${sample_id}.sam\n\n")
	elif endedness == "PE":
		a.write("bwa mem -M -t 10 " + bwa_dir + refgenome_prefix + " ${fq_r1} ${fq_r2} 2> " + \
				bwa_dir + prefix + "_${sample_id}_bwa-mem.out > " + \
				samtools_dir + prefix + "_${sample_id}.sam\n\n")

	a.write("samtools view -bS -F 4 " + samtools_dir + prefix + "_${sample_id}.sam > " + \
		bamtools_dir + prefix + "_${sample_id}.bam\n")
	a.write("rm " + samtools_dir + prefix + "_${sample_id}.sam\n\n")

	if endedness == "SE":
		a.write("samtools view -h -q 15 " + bamtools_dir + prefix + "_${sample_id}.bam" + \
			" | samtools view -buS - | samtools sort -o " + \
			bamtools_dir + prefix + "_${sample_id}_sorted.bam\n")
	elif endedness == "PE":
		a.write("samtools view -h " + bamtools_dir + prefix + "_${sample_id}.bam" + \
			" | samtools view -buS - | samtools sort -o " + \
			bamtools_dir + prefix + "_${sample_id}_sorted.bam\n")
	a.write("rm " + bamtools_dir + prefix + "_${sample_id}.bam\n\n")

	a.write("java -jar $PICARD MarkDuplicates" + \
			" I=" + bamtools_dir + prefix + "_${sample_id}_sorted.bam" + \
			" O=" + bamtools_dir + prefix + "_${sample_id}_sorted_dedup.bam" + \
			" M=" + bamtools_dir + prefix + "_${sample_id}_dups.log" + \
			" VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true\n")
	a.write("rm " + bamtools_dir + prefix + "_${sample_id}_sorted.bam\n\n")

	if endedness == "SE":
		a.write("samtools depth -aa " + bamtools_dir + prefix + "_${sample_id}_sorted_dedup.bam" + \
				" | cut -f 3 | gzip > " + bamtools_dir + prefix + "_${sample_id}.depth.gz")
		file_suffix = "_minq15_sorted_dedup.bam"
	elif endedness == "PE":
		a.write("bam clipOverlap --in " + bamtools_dir + prefix + "_${sample_id}_sorted_dedup.bam" + \
				" --out " + bamtools_dir + prefix + "_${sample_id}_sorted_dedup_clipped.bam --stats\n")
		a.write("rm " + bamtools_dir + prefix + "_${sample_id}_sorted_dedup.bam\n\n")

		a.write("samtools depth -aa " + bamtools_dir + prefix + "_${sample_id}_sorted_dedup_clipped.bam" + \
				" | cut -f 3 | gzip > " + bamtools_dir + prefix + "_${sample_id}.depth.gz")
		file_suffix = "_sorted_dedup_clipped.bam"

#Check depth with an array
depth_array_input = scripts_dir + prefix + "_depthsARRAY_input.txt"
iterator = 1
with open(depth_array_input, 'w') as i:
	for sample in trimmed_fastqs.keys():
		i.write(str(iterator) + ":" + bamtools_dir + prefix + "_" + sample + ".depth.gz\n")
		iterator += 1

#Write the depth calculations job array script
depth_array_script = scripts_dir + prefix + "_depthsARRAY.sh"
with open(depth_array_script, 'w') as d:
	d.write("#!/bin/bash\n\n")
	d.write("#SBATCH --job-name=depth\n")
	d.write("#SBATCH --cpus-per-task=5\n")
	d.write("#SBATCH --output=" + jobsout_dir + prefix + "depths_%A-%a.out" + "\n")
	d.write("#SBATCH --time=7-00:00:00\n")
	d.write("#SBATCH --array=1-" + str(iterator) + "%32\n\n")
	d.write("JOBS_FILE=" + depth_array_input + "\n")
	d.write("IDS=$(cat ${JOBS_FILE})\n\n")
	d.write("for sample_line in ${IDS}\n")
	d.write("do\n")
	d.write("""\tjob_index=$(echo ${sample_line} | awk -F ":" '{print $1}')\n""")
	d.write("""\tdepth_file=$(echo ${sample_line} | awk -F ":" '{print $2}')\n""")
	d.write("\tif [[ ${SLURM_ARRAY_TASK_ID} == ${job_index} ]]; then\n")
	d.write("\t\tbreak\n")
	d.write("\tfi\n")
	d.write("done\n\n")
	d.write("touch " + bamtools_dir + prefix + "_depths.csv\n")
	d.write("mean_cov_ind.py -i ${depth_file} -o " + bamtools_dir + prefix + "_depths.csv\n")

#Update checkpoint file with .bam filenames resulting from alignment
bamslist_file = working_dir + prefix + "_bamslist.txt"
with open(bamslist_file, 'w') as bams:
	for sample in trimmed_fastqs.keys():
		bams.write(bamtools_dir + prefix + "_" + sample + file_suffix + "\n")
with open(args.ckpt_file, 'a') as ckpt:
	ckpt.write("bamsLIST-all\t" + bamslist_file + "\n")

print("Step 3 has finished successfully! You will find two new scripts in ./scripts/: " + \
	align_array_script + " and " + depth_array_script + ".")
print("These scripts must be run in the order given above and both must finish before moving on to step 4.")
print("After " + depth_array_script + " has run, download the resulting file: " + prefix + "_depths.csv " + \
	"and generate a barchart of mean depths by individual. This will help you determine " + \
	"whether any individuals fall substantially below the average depth (usually, we use a cutoff of 1x).")
print("If you identify samples with coverage that is 'too low', add the sample id(s) to a new file, " + \
	"referred to as the 'blacklist' of individuals to be excluded from genotype likelhiood calculations " +\
	"and the final data sets.")
print("After generating this blacklist, you can continue to step 4 to write scripts for generating the final data sets.")
print("Remember to pass the checkpoint (.ckpt) file with the '-p' flag AND the blacklsit file with the '-b' flag.")
