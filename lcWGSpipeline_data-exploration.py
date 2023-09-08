#!/usr/bin/python3
import argparse
import os
from collections import OrderedDict

#Read in config file for the run
parser = argparse.ArgumentParser()
parser.add_argument('--config', '-c', help = 'A config file specifying the parameters for the run. This is not the pipeline config file!')
args = parser.parse_args()

working_dir = None
ref_genome = None
de_prefix = None
endedness = None
foldedness = None
fold_flag = None
minq = None
minmapq = None
min_depth = None
max_depth = None
trans_check = None
trans_options = None
bams_file = None
email = None

#Parse the config file to determine what's needed
with open(args.config, 'r') as c:
	for raw_line in c:
		param = raw_line.rstrip().split('\t')
		if param[0] == "working_directory":
			working_dir = param[1]
		elif param[0] == "reference_genome":
			ref_genome = param[1]
		elif param[0] == "data_exploration_prefix":
			de_prefix = param[1]
		elif param[0] == "endedness":
			endedness = param[1]
		elif param[0] == "foldedness":
			foldedness = param[1]
			if foldedness == "folded":
				fold_flag = "-fold 1 "
			else:
				fold_flag = ""
		elif param[0] == "minQ_vals":
			minq = param[1].split(',')
		elif param[0] == "minMapQ_vals":
			minmapq = param[1].split(',')
		elif param[0] == "min_depth":
			min_depth = param[1]
		elif param[0] == "max_depth":
			max_depth = param[1]
		elif param[0] == "check_transitions":
			if param[1] == "yes":
				trans_check = True
				trans_options = ["0", "1"]
			else:
				trans_check = False
				trans_options = ["0"]
		elif param[0] == "bamslist_file":
			bams_file = param[1]
		elif param[0] == "email":
			email = param[1]

#run a check: how many datasets will be generated if this moves forward?
#If this bit works, I'm going to be pretty excited...
combinatorics = len(minq) * len(minmapq) * len(trans_options)
continue_tf = None
continue_yn = None
extreme_continue_yn = None #allows me to throw a different user-check-in statement.

if combinatorics <= 4:
	continue_tf = True
elif 4 < combinatorics <= 8:
	while continue_yn not in ('yes', 'y', 'Yes', 'Y', 'no', 'n', 'No', 'N'):
		continue_yn = input("Testing all combinations of minQ, minMapQ, and including/excluding transitions " + \
			"will result in scripts to generate {0} datasets.\n".format(str(combinatorics)) + \
			"Do you wish to continue? (y/n)\n")
	if continue_yn in ('yes', 'y', 'Yes', 'Y'):
		continue_tf = True
	elif continue_yn in ('no', 'n', 'No', 'N'):
		continue_tf = False
else:
	while extreme_continue_yn not in ('yes', 'y', 'Yes', 'Y', 'no', 'n', 'No', 'N'):
		extreme_continue_yn = input("Testing all combinations of minQ, minMapQ, and including/excluding transitions " + \
			"will result in scripts to generate {0} datasets.\n".format(str(combinatorics)) + \
			"This is a lot. You may continue, but I recommend running these scripts in series, rather than parallel.\n" + \
			"Do you wish to continue? (y/n)\n")
	if extreme_continue_yn in ('yes', 'y', 'Yes', 'Y'):
		continue_tf = True
	elif extreme_continue_yn in ('no', 'n', 'No', 'N'):
		continue_tf = False

if continue_tf:
	#setup the environment a little bit
	scripts_dir = working_dir + "scripts/"
	if os.path.isdir(scripts_dir) is not True:
		os.mkdir(scripts_dir)
	jobsout_dir = working_dir + "jobsout/"
	if os.path.isdir(jobsout_dir) is not True:
		os.mkdir(jobsout_dir)
	results_dir = working_dir + "results/"
	if os.path.isdir(results_dir) is not True:
		os.mkdir(results_dir)

	#generate the parameter combinations for heterozygosity estimation
	parameterizations = OrderedDict()
	for mq in minq:
		for mmq in minmapq:
			for trans in trans_options:
				# this is a little counterintuitive and I might rework the transitions bit, -noTrans takes 0 to ALLOW transitions and 1 to DISALLOW transitions
				run_suffix = "mindp" + min_depth + "_maxdp" + max_depth + "_minq" + mq + "_minmapq" + mmq + "_notransitions" + trans + "_" + foldedness
				parameterizations[run_suffix] = fold_flag + "-minQ " + mq + " -minMapQ " + mmq + " -noTrans " + trans

	#write the input file for all arrays
	de_array_input = scripts_dir + de_prefix + "_heterozygosityARRAY_input.txt"
	iterator = 0
	with open(bams_file, 'r') as i:
		with open(de_array_input, 'w') as o:
			for raw_bamfilename in i:
				iterator += 1 #this way, the iterator is the array length, without needing to subtract 1 later
				o.write(str(iterator) + ":" + raw_bamfilename)

	for angsd_run_title, angsd_run_specifics in parameterizations.items():
		angsd_run_script = scripts_dir + de_prefix + "_" + angsd_run_title + "_heterozygosityARRAY.sh"
		with open(angsd_run_script, 'w') as s:
			s.write("#!/bin/bash\n\n")
			s.write("#SBATCH --cpus-per-task=8\n")
			s.write("#SBATCH --time=1-00:00:00\n") #re-evaluate the appropriateness of this value after getting some IRL
			s.write("#SBATCH --job-name=de\n")
			s.write("#SBATCH --output=" + jobsout_dir + de_prefix + "_" + angsd_run_title + "_%A-%a.out\n")
			s.write("#SBATCH --mail-type=FAIL\n")
			s.write("#SBATCH --mail-user=" + email + "\n")
			s.write("#SBATCH --array=1-" + str(iterator) + "%24\n\n")
			s.write("module unload bio/angsd/0.933\n")
			s.write("module load bio/angsd/0.933\n\n")

			s.write("JOBS_FILE=" + de_array_input + "\n")
			s.write("IDS=$(cat ${JOBS_FILE})\n\n")

			s.write("for sample_line in ${IDS}\n")
			s.write("do\n")
			s.write("""\tjob_index=$(echo ${sample_line} | awk -F ":" '{print $1}')\n""")
			s.write("""\tbamfile=$(echo ${sample_line} | awk -F ":" '{print $2}')\n""")
			s.write("\tif [[ ${SLURM_ARRAY_TASK_ID} == ${job_index} ]]; then\n")
			s.write("\t\tbreak\n")
			s.write("\tfi\n")
			s.write("done\n\n")

			s.write("sample_id=$(echo $bamfile | sed 's!^.*/!!')\n")
			s.write("sample_id=$(echo $sample_id | awk -F _ '{print $2}')\n\n")

#get SAF file
			s.write("angsd " + \
				"-i ${bamfile} " + \
				"-anc " + ref_genome + " " + \
				"-ref " + ref_genome + " " + \
				"-out " + results_dir + de_prefix + "_${sample_id}_" + angsd_run_title + " " + \
				"-nThreads 8 " + \
				"-doSaf 1 " + \
				"-GL 1 " + \
				"-doCounts 1 " + \
				"-setMinDepth " + min_depth + " " + \
				"-setMaxDepth " + max_depth + " " + \
				"-C 50 " + \
				"-remove_bads 1 " + \
				"-uniqueOnly 1 " + \
				"-trim 0 " + \
				angsd_run_specifics)
			if endedness == "PE":
				s.write(" -only_proper_pairs 1\n\n")
			else:
				s.write("\n\n")

#estimate SFS
#tole is tolerance for breaking EM (When the difference in successive likelihood values in the EM algorithm gets below this value the optimization will stop)
	#0.0000001 comes from Nicolas' code
			s.write("realSFS " + \
				results_dir + de_prefix + "_${sample_id}_" + angsd_run_title + ".saf.idx " + \
				"-tole 0.0000001 " + \
				"-P 8 " + \
				"-seed 37 " + \
				"> " + results_dir + de_prefix + "_${sample_id}_" + angsd_run_title + ".ml\n\n")

			s.write("touch " + results_dir + de_prefix + "_" + angsd_run_title + ".het\n")
			s.write("echo -ne ${sample_id}'\t' >> " + results_dir + de_prefix + "_" + angsd_run_title + ".het\n")
			if foldedness == "folded":
				s.write("awk '{ val1=$1 ; val2=$2 ; het_val=val2/(val1+val2); print het_val; }' ")
			else:
				s.write("awk '{ val1=$1 ; val2=$2 ; val3=$3 ; het_val=val2/(val1+val2+val3); print het_val; }' ")
			s.write(results_dir + de_prefix + "_${sample_id}_" + angsd_run_title + ".ml >> " + \
				results_dir + de_prefix + "_" + angsd_run_title + ".het\n")

	print("Scripts have been written for data exploration.")
	print("The execution of each script will result in a file of heterozygosity values for every individual included in the user-specified bamfile.")
else:
	print("You have opted to stop. No scripts have been written.")