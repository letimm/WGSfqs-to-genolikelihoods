#!/usr/bin/python3
import argparse
import subprocess

#Read in config file for the run
parser = argparse.ArgumentParser()
parser.add_argument('--ckpt_file', '-p', help = 'Please provide the checkpoint file created in step 0.')
args = parser.parse_args()

#Initialize run config variables with some default values
scripts_dir = None
bwa_dir = None
samtools_dir = None
bamtools_dir = None
n_nodes = None
prefix = None
endedness = None
samslist = []
sample_ids = []
bamslist = []

#Parse the config file to determine what's needed (if user wants bwa, no need for picard or bowtie2)
with open(args.ckpt_file, 'r') as last_step_ckpt:
	for raw_ckpt_line in last_step_ckpt:
		ckpt_line = raw_ckpt_line.rstrip()
		ckpt_setting = ckpt_line.split('\t')
		if ckpt_setting[0] == "scriptsDIR":
			scripts_dir = ckpt_setting[1]
		elif ckpt_setting[0] == "bwaDIR":
			bwa_dir = ckpt_setting[1]
		elif ckpt_setting[0] == "samtoolsDIR":
			samtools_dir = ckpt_setting[1]
		elif ckpt_setting[0] == "bamtoolsDIR":
			bamtools_dir = ckpt_setting[1]
		elif ckpt_setting[0] == "distNODES":
			n_nodes = float(ckpt_setting[1])
		elif ckpt_setting[0] == "prefix":
			prefix = ckpt_setting[1]
		elif ckpt_setting[0] == "ENDEDNESS":
			endedness = ckpt_setting[1]
		elif ckpt_setting[0] == "SAM":
			samslist.append(ckpt_setting[1])

#Take each .sam file through SAMTOOLS to generate a quality-filtered .bam file
for samfile in samslist:
	samfile_name_as_list = samfile.split(".")
	sample_id = ".".join(samfile_name_as_list[:-1])
	first_bam = sample_id + ".bam"
	filtered_bam = sample_id + "_minq20.bam"
	sorted_bam = sample_id + "_minq20_sorted.bam"
	dedup_bam = sample_id + "_minq20_sorted_dedup.bam"
	clipped_bam = sample_id + "_minq20_sorted_dedup_clipped.bam"
	duplicates_log = sample_id + "_dups.log"
	script = scripts_dir + sample_id + "_sambamSLURM.sh"
	with open(script, 'w') as s:
		s.write("#!/bin/bash\n\n")
		s.write("#SBATCH --nodes=" + str(int(n_nodes)) + "\n")
		s.write("#SBATCH --ntasks=" + str(int(n_nodes)) + "\n")
		s.write("#SBATCH --job-name=sam_" + sample_id + "\n")
		s.write("#SBATCH --output=" + samtools_dir + sample_id + "_sambam.out\n\n")
		s.write("module unload bio/samtools bio/picard bio/bamtools\n")
		s.write("module load bio/samtools bio/picard bio/bamtools\n\n")
		s.write("samtools fixmate -O bam " + samtools_dir + samfile + " " + bamtools_dir + first_bam + "\n")
		s.write("samtools view -b -F 4 " + bamtools_dir + first_bam + " > " + bamtools_dir + filtered_bam + "\n")
		s.write("samtools view -h -q 20 " + bamtools_dir + filtered_bam + \
			" | samtools view -buS - | samtools sort -o " + bamtools_dir + sorted_bam + "\n")
		s.write("java -jar $PICARD MarkDuplicates \
			I=" + bamtools_dir + sorted_bam + " O=" + bamtools_dir + dedup_bam + \
			" M=" + bamtools_dir + duplicates_log + \
			" VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true\n")
		if endedness == "SE":
			bamslist.append(dedup_bam)
			s.write("samtools depth -aa " + bamtools_dir + dedup_bam + \
				" | cut -f 3 | gzip > " + bamtools_dir + sample_id + ".depth.gz\n")
		elif endedness == "PE":
			bamslist.append(clipped_bam)
			s.write("bamtools clipOverlap --in " + bamtools_dir + dedup_bam + \
				" --out " + bamtools_dir + clipped_bam + "--stats\n")
			s.write("samtools depth -aa " + bamtools_dir + clipped_bam + \
				" | cut -f 3 | gzip > " + bamtools_dir + clipped_bam + ".depth.gz")
	subprocess.call(["sbatch", script])

#Update checkpoint file with .bam filenames resulting from SAMTOOLS run
bams_list_file = prefix + "_bams_list.txt"
with open(bams_list_file, 'w') as bams:
	for bamfile in bamslist:
		bams.write(bamtools_dir + bamfile + "\n")

with open(args.ckpt_file, 'a') as ckpt:
	ckpt.write("bamsLIST\t" + bams_list_file +"\n")