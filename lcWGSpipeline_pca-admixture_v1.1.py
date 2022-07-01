#!/usr/bin/python3
import argparse
import gzip
import math
import os
from collections import OrderedDict

#Read in config file for the run
parser = argparse.ArgumentParser()
parser.add_argument('--ckpt_file', '-p', help = 'Please provide the checkpoint file created in step 0.')
parser.add_argument('--k_val_max', '-k', default = 10, help = 'Please provide the maximum k value you would like to test.')
args = parser.parse_args()

#Initialize run config variables with some default values
working_dir = None
scripts_dir = None
jobsout_dir = None
gls_dir = None
prefix = None
polymorphic_beagle_files = []
wgp_beagle = None

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
		elif ckpt_setting[0] == "glsDIR":
			gls_dir = ckpt_setting[1]
		elif ckpt_setting[0] == "prefix":
			prefix = ckpt_setting[1]
		elif ckpt_setting[0] == "glsFILES-polymorphic":
			polymorphic_beagle_files = ckpt_setting[1].split(",")

pca_dir = working_dir + "pca/"
admix_dir = working_dir + "admixture/"
if os.path.isdir(pca_dir) is not True:
	os.mkdir(pca_dir)
if os.path.isdir(admix_dir) is not True:
	os.mkdir(admix_dir)

pca_array_input = scripts_dir + prefix + "_pcangsdARRAY_input.txt"
with open(pca_array_input, 'w') as i:
	iterator = 1
	for beagle in polymorphic_beagle_files:
		i.write(str(iterator) + ":" + beagle + "\n")
		iterator += 1

pcangsd_script = scripts_dir + prefix + "_pcangsdARRAY.sh"
with open(pcangsd_script, 'w') as pca:
	pca.write("#!/bin/bash\n\n")
	pca.write("#SBATCH --cpus-per-task=10\n")
	pca.write("#SBATCH --time=7-00:00:00\n")
	pca.write("#SBATCH --job-name=pca_" + prefix + "\n")
	pca.write("#SBATCH --output=" + jobsout_dir + prefix + "_pcangsd_%A-%a.out\n")
	pca.write("#SBATCH --array=1-" + str(len(polymorphic_beagle_files)) + "%24\n\n")
	pca.write("module unload bio/pcangsd/0.99\n")
	pca.write("module load bio/pcangsd/0.99\n")
	pca.write("source /opt/bioinformatics/venv/pcangsd-0.99/bin/activate\n\n")

	pca.write("JOBS_FILE=" + pca_array_input + "\n")
	pca.write("IDS=$(cat ${JOBS_FILE})\n\n")
	pca.write("for sample_line in ${IDS}\n")
	pca.write("do\n")
	pca.write("""\tjob_index=$(echo ${sample_line} | awk -F ":" '{print $1}')\n""")
	pca.write("""\tbeagle_file=$(echo ${sample_line} | awk -F ":" '{print $2}')\n""")
	pca.write("\tif [[ ${SLURM_ARRAY_TASK_ID} == ${job_index} ]]; then\n")
	pca.write("\t\tbreak\n")
	pca.write("\tfi\n")
	pca.write("done\n\n")

	pca.write("chrom=$(echo $beagle_file | sed 's!^.*/!!')\n")
	pca.write("chrom=${chrom%.beagle.gz}\n\n")

	pca.write("pcangsd.py" + \
		" -threads 10" + \
		" -beagle ${beagle_file}" + \
		" -o " + pca_dir + "${chrom}" + \
		" -sites_save" + \
		" -pcadapt")

#Concatenate beagle files to make the gls file for the whole genome
all_chrs_poly_beagle = gls_dir + prefix + "_wholegenome_polymorphic.beagle"
beagle_concat_script = scripts_dir + prefix + "_concatenate_beagles.sh"
with open(beagle_concat_script, 'w') as cs:
	cs.write("#!/bin/bash\n\n")
	cs.write("#SBATCH --cpus-per-task=10\n")
	cs.write("#SBATCH --time=7-00:00:00\n")
	cs.write("#SBATCH --job-name=" + prefix + "_concat-beagles\n")
	cs.write("#SBATCH --output=" + jobsout_dir + prefix + "_concatenate-beagles_%A.out\n\n")
	cs.write("zcat " + polymorphic_beagle_files[0] + " " + \
		"| head -n 1 > " + \
		all_chrs_poly_beagle + "; " + \
		"for i in " + " ".join(polymorphic_beagle_files) + "; " + \
		"do zcat $i " + \
		"| tail -n +2 -q >> " + \
		all_chrs_poly_beagle + "; " + \
		"done\n")
	cs.write("gzip " + all_chrs_poly_beagle)

whole_genome_pca_script = scripts_dir + prefix + "_wholegenome_pcangsd.sh"
with open(whole_genome_pca_script, 'w') as wgp:
	wgp.write("#!/bin/bash\n\n")
	wgp.write("#SBATCH --cpus-per-task=10\n")
	wgp.write("#SBATCH --time=7-00:00:00\n")
	wgp.write("#SBATCH --job-name=" + prefix + "_wgp-pca-admix\n")
	wgp.write("#SBATCH --output=" + jobsout_dir + prefix + "_wholegenome_polymorphic_%A.out\n\n")
	wgp.write("module unload bio/pcangsd/0.99\n")
	wgp.write("module load bio/pcangsd/0.99\n")
	wgp.write("source /opt/bioinformatics/venv/pcangsd-0.99/bin/activate\n\n")

	wgp.write("pcangsd.py " + \
		"-threads 10 " + \
		"-beagle " + all_chrs_poly_beagle + ".gz " + \
		"-o " + pca_dir + prefix + "_wholegenome-polymorphic " + \
		"-sites_save " + \
		"-pcadapt\n\n")

#Make an array script that tests all K values 1-max_k
whole_genome_admix_script = scripts_dir + prefix + "_wholegenome_admixARRAY.sh"
with open(whole_genome_admix_script, 'w') as wga:
	wga.write("#!/bin/bash\n\n")
	wga.write("#SBATCH --cpus-per-task=10\n")
	wga.write("#SBATCH --time=7-00:00:00\n")
	wga.write("#SBATCH --job-name=" + prefix + "_wgp-pca-admix\n")
	wga.write("#SBATCH --output=" + jobsout_dir + prefix + "_wholegenome_polymorphic_%A-%a.out\n\n")
	wga.write("#SBATCH --array=1-" + args.k_val_max + "%12\n\n")

	wga.write("module unload bio/ngsadmix\n")
	wga.write("module load bio/ngsadmix\n\n")

	wga.write("for k_val in {1.." + args.k_val_max + "}\n")
	wga.write("do\n")
	wga.write("\tif [[ ${SLURM_ARRAY_TASK_ID} == ${k_val} ]]; then\n")
	wga.write("\t\tbreak\n")
	wga.write("\tfi\n")
	wga.write("done\n\n")

	wga.write("NGSadmix " + \
		"-likes " + all_chrs_poly_beagle + ".gz " + \
		"-K ${k_val} " + \
		"-outfiles " + admix_dir + prefix + "_wholegenome-polymorphic_k${k_val} " + \
		"-P 10 " + \
		"-minMaf 0")

with open(args.ckpt_file, 'a') as ckpt:
	ckpt.write("all-chrs_glsFILE-polymorphic\t" + all_chrs_poly_beagle + ".gz\n")
	ckpt.write("pcaDIR\t" + pca_dir + "\n")
	ckpt.write("admixtureDIR\t" + admix_dir + "\n")

print("Four scripts have been generated: " + beagle_concat_script + ", " + pcangsd_script + ", " + whole_genome_pca_script + ", and " + whole_genome_admix_script + ".")
print(beagle_concat_script + " must run before " + whole_genome_pca_script + " and " + whole_genome_admix_script + ".")
print(pcangsd_script + " will run pca for each chromosome.")
print(whole_genome_pca_script + " will run pca for the whole genome polymorphic data.")
print(whole_genome_admix_script + " will launch an array of admixture jobs, testing every value of K between 1 and the max k value entered by the user.")
