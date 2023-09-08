#!/usr/bin/python3
import argparse
from collections import OrderedDict

#Read in config file for the run
parser = argparse.ArgumentParser()
parser.add_argument('--ckpt_file', '-p', help = 'Please provide the checkpoint file created in step 0.')
parser.add_argument('--downsample_prefix', '-d', help = 'Option to provide a prefix for the dataset/results including downsampled individuals (if you intend to test across different aspects of your data (batch, region, population), I suggest setting this, rather than running with the defualt. Default: PREFIXds', required = False)
parser.add_argument('--array_input_file', '-a', help = 'Please provide the array input file generated with depths_check.R')
args = parser.parse_args()

#Initialize run config variables with some default values
working_dir = None
scripts_dir = None
jobsout_dir = None
prefix = None
email = None
bamtools_dir = None
bamslist_file = None

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
		elif ckpt_setting[0] == "prefix":
			prefix = ckpt_setting[1]
		elif ckpt_setting[0] == "email":
			email = ckpt_setting[1]
		elif ckpt_setting[0] == "bamtoolsDIR":
			bamtools_dir = ckpt_setting[1]
		elif ckpt_setting[0] == "bamsLIST-all":
			full_bams_list = ckpt_setting[1]

# First, get a list of downsampled files
downsampled_individuals = []
with open(args.array_input_file, 'r') as di:
	#just doing all this in one line, but can split it out if too hard to look at
	for line in di:
		downsampled_individuals.append(line.rstrip().split(':')[1])

downsample_script = scripts_dir + prefix + "_downsample-bamsARRAY.sh"
with open(downsample_script, 'w') as ds:
	ds.write("#!/bin/bash\n\n")
	ds.write("#SBATCH --cpus-per-task=5\n")
	ds.write("#SBATCH --time=0-12:00:00\n")
	ds.write("#SBATCH --job-name=ds_" + prefix + "\n")
	ds.write("#SBATCH --output=" + jobsout_dir + prefix + "_downsample-bams_%A-%a.out\n")
	ds.write("#SBATCH --mail-type=FAIL\n")
	ds.write("#SBATCH --mail-user=" + email + "\n")
	ds.write("#SBATCH --array=1-" + str(len(downsampled_individuals)) + "%48\n\n")
	ds.write("module unload bio/samtools/1.11\n")
	ds.write("module load bio/samtools/1.11\n\n")

	ds.write("JOBS_FILE=" + args.array_input_file + "\n")
	ds.write("IDS=$(cat ${JOBS_FILE})\n\n")
	ds.write("for sample_line in ${IDS}\n")
	ds.write("do\n")
	ds.write("""\tjob_index=$(echo ${sample_line} | awk -F ":" '{print $1}')\n""")
	ds.write("""\tsample_id=$(echo ${sample_line} | awk -F ":" '{print $2}')\n""")
	ds.write("""\tdownsample_value=$(echo ${sample_line} | awk -F ":" '{print $3}')\n""")
	ds.write("\tif [[ ${SLURM_ARRAY_TASK_ID} == ${job_index} ]]; then\n")
	ds.write("\t\tbreak\n")
	ds.write("\tfi\n")
	ds.write("done\n\n")

	ds.write("samtools view " + \
		"-bo " + bamtools_dir + prefix + "_${sample_id}_sorted_dedup_clipped_downsampled.bam " + \
		"-s ${downsample_value} " + \
		bamtools_dir + prefix + "_${sample_id}_sorted_dedup_clipped.bam\n")
	ds.write("samtools depth " + \
		"-aa " + bamtools_dir + prefix + "_${sample_id}_sorted_dedup_clipped_downsampled.bam " + \
		"| cut -f 3 | " + \
		"gzip > " + bamtools_dir + prefix + "_${sample_id}_sorted_dedup_clipped_downsampled.depth.gz\n")
	ds.write("samtools index " + \
		bamtools_dir + prefix + "_${sample_id}_sorted_dedup_clipped_downsampled.bam")

#Update the prefix for the dataset containing downsampled inds
if args.downsample_prefix == None:
	downsample_prefix = prefix + "-ds"
else:
	downsample_prefix = prefix + "-" +  args.downsample_prefix

# write out the new bamslist
ds_bamslist_file = working_dir + downsample_prefix + "_bamslist.txt"
with open(full_bams_list, 'r') as bi:
	with open(ds_bamslist_file, 'w') as bo:
		for raw_line in bi:
			sample_bam = raw_line.rstrip()
			sample_id = sample_bam.split('/')[-1].split('_')[1]
			if sample_id in downsampled_individuals:
				bo.write(sample_bam.split('.')[0] + "_downsampled.bam\n")
			else:
				bo.write(raw_line)

#Spawn new checkpoints
qc_list = ["15", "25", "35"]
ckpt_files = []
for qc_val in qc_list:
	dsq_ckpt_filename = working_dir + downsample_prefix + "-q" + qc_val + ".ckpt"
	ckpt_files.append(dsq_ckpt_filename)
	with open(args.ckpt_file, 'r') as ci:
		with open(dsq_ckpt_filename, 'w') as co:
			for line in ci:
				if line.rstrip().split('\t')[0] == "prefix":
					co.write("prefix\t" + downsample_prefix + "-q" + qc_val + "\n")
				elif line.rstrip().split('\t')[0] == "bamsLIST-all":
					co.write("bamsLIST-all\t" + ds_bamslist_file + "\n")
				else:
					co.write(line)
			co.write("qFILTER\t" + qc_val + "\n")

#Update checkpoint file with the names of spawned ckpts
with open(args.ckpt_file, 'a') as ckpt:
	ckpt.write("branchedCKPTs\t" + ",".join(ckpt_files) + "\n")

print("This optional step has run successfully! A new script in ./scripts/:\n" + \
	downsample_script + " will downsample all individuals identified by depth analysis.\n")
print("Additionally, three new ckpt files have been generated: " + \
	"\n".join(ckpt_files) + "\n"
	"Call lcWGSpipeline_step4 on each of these scripts, as desired for data exploration.")
