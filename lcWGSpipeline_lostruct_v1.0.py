#!/usr/bin/python3
import argparse
import gzip
import math
import os
from collections import OrderedDict

#Read in config file for the run
parser = argparse.ArgumentParser()
parser.add_argument('--ckpt_file', '-p', help = 'Please provide the checkpoint file created in step 0.')
parser.add_argument('--window_size', '-w', default = "100000", help = 'Provide the desired window size.')
parser.add_argument('--num_pcs', '-n', default = "3", help = 'Provide the number of principal components to retain for each window.')
parser.add_argument('--rscript_loc', '-r', help = 'Provide the full path to the R script that converts cov to eigenvalues (if you grab the R script from the same repo you found this, you want "cov_to_eigen.R").')
args = parser.parse_args()

#Initialize run config variables with some default values
working_dir = None
scripts_dir = None
jobsout_dir = None
gls_dir = None
prefix = None
beagles = []

#Parse the config file
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
		elif ckpt_setting[0] == "glsFILES-polymorphic":
			beagles = ckpt_setting[1].split(",")

#Make lostruct directory
lostruct_dir = working_dir + "lostruct/"
if os.path.isdir(lostruct_dir) is not True:
	os.mkdir(lostruct_dir)

#Grab beagle header
beagle_header = None
with gzip.open(beagles[0], 'rt') as abf:
	beagle_header = abf.readline()

#Make an array script
eigen_array_input = scripts_dir + prefix + "_eigenARRAY_input.txt"
with open(eigen_array_input, 'w') as e:
	for i, beagle_file in enumerate(beagles):
		iterator = i + 1
		beagle_prefix = lostruct_dir + beagle_file.split("/")[-1].rstrip(".beagle.gz")
		chromosome = beagle_file.split("_")[-2]
		e.write(str(iterator) + ":" + beagle_file + ":" + beagle_prefix + ":" + chromosome + "\n")

eigen_script = scripts_dir + prefix + "_eigenARRAY.sh"
with open(eigen_script, 'w') as eig:
	eig.write("#!/bin/bash\n\n")
	eig.write("#SBATCH --cpus-per-task=5\n")
	eig.write("#SBATCH --time=7-00:00:00\n")
	eig.write("#SBATCH --job-name=eigen_" + prefix + "\n")
	eig.write("#SBATCH --output=" + jobsout_dir + prefix + "_eigen_%A-%a.out\n")
	eig.write("#SBATCH --array=1-" + str(len(beagles)) + "%48\n\n")
	eig.write("module unload bio/pcangsd/0.99 R/4.0.3\n")
	eig.write("module load bio/pcangsd/0.99 R/4.0.3\n")
	eig.write("source /opt/bioinformatics/venv/pcangsd-0.99/bin/activate\n\n")

	eig.write("JOBS_FILE=" + eigen_array_input + "\n")
	eig.write("IDS=$(cat ${JOBS_FILE})\n\n")
	eig.write("for sample_line in ${IDS}\n")
	eig.write("do\n")
	eig.write("""\tjob_index=$(echo ${sample_line} | awk -F ":" '{print $1}')\n""")
	eig.write("""\tbeagle=$(echo ${sample_line} | awk -F ":" '{print $2}')\n""")
	eig.write("""\tsplit_beagle_prefix=$(echo ${sample_line} | awk -F ":" '{print $3}')\n""")
	eig.write("""\tchrom=$(echo ${sample_line} | awk -F ":" '{print $4}')\n""")
	eig.write("\tif [[ ${SLURM_ARRAY_TASK_ID} == ${job_index} ]]; then\n")
	eig.write("\t\tbreak\n")
	eig.write("\tfi\n")
	eig.write("done\n\n")

	eig.write("zcat ${beagle} | tail -n +2 > ${beagle%.gz}\n")
	eig.write("split -l " + args.window_size + " ${beagle%.gz} ${split_beagle_prefix}_part-\n")
	eig.write("rm ${beagle%.gz}\n")
	eig.write("touch ${split_beagle_prefix}_puppy-eigens.txt\n")
	eig.write("for i in ${split_beagle_prefix}_part-*\n")
	eig.write("do\n")
	eig.write("\tsed -i '1s/^/" + beagle_header.replace("\n","\\n") + "/' $i\n")
	eig.write("\tgzip $i\n")
	eig.write("\tpcangsd.py" + \
		" -threads 5" + \
		" -beagle $i.gz" + \
		" -o $i" + \
		" -sites_save" + \
		" -pcadapt\n")
	eig.write("\tRscript --vanilla " + \
		args.rscript_loc + " " + \
		"$i.cov " + \
		args.num_pcs + " " + \
		args.window_size + " " + \
		"$i.gz " + \
		"${chrom} " + \
		lostruct_dir + "\n")
#	eig.write("\t" append the rscript outfile to puppy-eigens.txt "\n")
	eig.write("done")


#Update checkpoint file with the lostruct directory and the eigen files
#with open(args.ckpt_file, 'a') as ckpt:
#	ckpt.write("lostructDIR\t" + lostruct_dir + "\n")
#	ckpt.write("eigenFILES\t" +  + "\n")
