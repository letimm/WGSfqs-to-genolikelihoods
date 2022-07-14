#!/usr/bin/python3
import argparse
import os
from collections import OrderedDict

#Read in config file for the run
parser = argparse.ArgumentParser()
parser.add_argument('--ckpt_file', '-p', help = 'Please provide the checkpoint file created in step 0.')
parser.add_argument('--rscript_loc', '-r', help = 'Provide the full path to the HMM R script you would like to run.')
args = parser.parse_args()

#Initialize run config variables with some default values
working_dir = None
scripts_dir = None
jobsout_dir = None
fst_dir = None
prefix = None
fsts = OrderedDict()

#Parse the ckpt file
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
		elif ckpt_setting[0] == "fstDIR":
			fst_dir = ckpt_setting[1]
		elif ckpt_setting[0] == "prefix":
			prefix = ckpt_setting[1]
		elif ckpt_setting[0] == "population-pair_fst-popFILES":
			fsts[ckpt_setting[1]] = ckpt_setting[2].split(",")

#make hmm directory
hmm_dir = working_dir + "hmm/"
if os.path.isdir(hmm_dir) is not True:
	os.mkdir(hmm_dir)

#write out a new file of fst values only (formatted for input to hmm)
hmm_infiles = OrderedDict()
hmm_detailed_infiles = OrderedDict()
hmm_outfiles = OrderedDict()
num_jobs = 0
for poppair, poppair_fsts in fsts.items():
	#this file records the fst values for all sites (just fsts - will be input for hmm in R)
	hmm_in = hmm_dir + prefix + "_wholegenome_" + poppair + "_fst-hmm-in.txt"
	hmm_infiles[poppair] = hmm_in
	#this file records the details for all sites (details!)
	hmm_detail_in = hmm_dir + prefix + "_wholegenome_" + poppair + "_fst-hmm-detail.txt"
	hmm_detailed_infiles[poppair] = hmm_detail_in
	#not writing anything out here, just recording the hmm outfile name; also counting how many files will be in the job array
	hmm_out = hmm_dir + prefix + "_wholegenome_" + poppair + "_fst-hmm-in_3state_HMMstates.txt"
	hmm_outfiles[poppair] = hmm_out
	num_jobs += 1

	with open(hmm_in, 'w') as h, open(hmm_detail_in, 'w') as d:
		h.write('"x"\n')
		d.write("region\tchr\tmidPos\tNsites\tfst\n")
		#retain fsts from all sites across the genome for the population pair
		for chr_fst in poppair_fsts:
			with open(chr_fst, 'r') as f:
				next(f)
				for raw_line in f:
					d.write(raw_line)
					fst_val = raw_line.rstrip().split("\t")[-1]
					h.write(fst_val + "\n")

hmm_array_input = scripts_dir + prefix + "_hmmARRAY_input.txt"
with open(hmm_array_input, 'w') as i:
	iterator = 1
	for hmm_file in hmm_infiles.values():
		i.write(str(iterator) + ":" + hmm_file + "\n")
		iterator += 1

hmm_script = scripts_dir + prefix + "_hmmARRAY.sh"
with open(hmm_script, 'w') as hmm:
	hmm.write("#!/bin/bash\n\n")
	hmm.write("#SBATCH --cpus-per-task=10\n")
	hmm.write("#SBATCH --time=7-00:00:00\n")
	hmm.write("#SBATCH --job-name=hmm_" + prefix + "\n")
	hmm.write("#SBATCH --output=" + jobsout_dir + prefix + "_hmm_%A-%a.out\n")
	hmm.write("#SBATCH --array=1-" + str(num_jobs) + "%12\n\n")
	hmm.write("module unload R/4.0.3\n")
	hmm.write("module load R/4.0.3\n\n")

	hmm.write("JOBS_FILE=" + hmm_array_input + "\n")
	hmm.write("IDS=$(cat ${JOBS_FILE})\n\n")
	hmm.write("for sample_line in ${IDS}\n")
	hmm.write("do\n")
	hmm.write("""\tjob_index=$(echo ${sample_line} | awk -F ":" '{print $1}')\n""")
	hmm.write("""\tfst_file=$(echo ${sample_line} | awk -F ":" '{print $2}')\n""")
	hmm.write("\tif [[ ${SLURM_ARRAY_TASK_ID} == ${job_index} ]]; then\n")
	hmm.write("\t\tbreak\n")
	hmm.write("\tfi\n")
	hmm.write("done\n\n")

	hmm.write("Rscript --vanilla " + args.rscript_loc + " ${fst_file} 10")

with open(args.ckpt_file, 'a') as ckpt:
	ckpt.write("hmmDIR\t" + hmm_dir + "\n")
	for two_pops, infile in hmm_detailed_infiles.items():
		ckpt.write("hmm-inFILES\t" + two_pops + "\t" + infile + "\n")
		ckpt.write("hmm-outFILES\t" + two_pops + "\t" + hmm_outfiles[two_pops] + "\n")

print("A new script has been generated to run HMM in R. Submitting " + hmm_script + " will classify each site as high (1), low (2), or background differentiation (3), based on the fst value\n")
print("Note: for this step to run correctly, you must have HMM_log10FST+1_3norm.R available for the script to access (in your path).\n")
