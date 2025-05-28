#!/usr/bin/python3

import argparse
import glob
import math
import os
import subprocess

#Read in config file for the run
parser = argparse.ArgumentParser()
parser.add_argument('--ckpt_file', '-p', help = 'Please provide the checkpoint file created in step 0.')
parser.add_argument('--alpha', '-a', type = str, help = 'Please provide significance threshold for the Chi-Square tests.', default = 0.05)
parser.add_argument('--parallelize_n', '-n', type = float, help = 'Please provide a number of parallel jobs to array (must be <999).', default = 100)
parser.add_argument('--ngsParalog_executable', '-x', help = 'Please provide the ngsParalog executable (include the full path, like: /home/ltimm/bin/ngsParalog/ngsParalog).')
parser.add_argument('--sigTest', '-s', help = 'Please provide the full path to ngsParalog_sigTest.R.')
args = parser.parse_args()

#Initialize run config variables with some default values
working_dir = None
scripts_dir = None
jobsout_dir = None
gls_dir = None
endedness = None
ref_genome = None
filtered_bamslist_file = None
filtered_n = 0
prefix = None
orig_prefix = None
email = None
sites_file = None
min_depth = None
max_depth = None
minq = None
min_ind = None
chrs_list = []

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
		elif ckpt_setting[0] == "usePREFIX":
			prefix = ckpt_setting[1]
		elif ckpt_setting[0] == "prefix":
			orig_prefix = ckpt_setting[1]
		elif ckpt_setting[0] == "email":
			email = ckpt_setting[1]
		elif ckpt_setting[0] == "ENDEDNESS":
			endedness = ckpt_setting[1]
		elif ckpt_setting[0] == "glsDIR":
			gls_dir = ckpt_setting[1]
		elif ckpt_setting[0] == "refgenomeFASTA":
			ref_genome = ckpt_setting[1]
		elif ckpt_setting[0] == "bamsLIST-filtered":
			filtered_bamslist_file = ckpt_setting[1]
		elif ckpt_setting[0] == "nIND-filtered":
			filtered_n = float(ckpt_setting[1])
		elif ckpt_setting[0] == "chrsLIST":
			chrs_list = ckpt_setting[1].split(",")
		elif ckpt_setting[0] == "wholegenomeSITES-polymorphic":
			sites_file = ckpt_setting[1]
#stuff for angsd
		elif ckpt_setting[0] == "depthMIN-factor":
			min_depth = float(ckpt_setting[1])
		elif ckpt_setting[0] == "depthMAX-factor":
			max_depth = float(ckpt_setting[1])
		elif ckpt_setting[0] == "minQ":
			minq = ckpt_setting[1]
		elif ckpt_setting[0] == "minIND":
			min_ind = ckpt_setting[1]

paralog_dir = working_dir + "paralog/"
if os.path.isdir(paralog_dir) is not True:
	os.mkdir(paralog_dir)

#given x, how many sites go to each job?
#how many sites in the file?
sites_ct = 0
with open(sites_file, "rb") as f:
    sites_ct = sum(1 for _ in f)

#how many lines per job (always round up)?
if args.parallelize_n > 999:
	real_n = 999
else:
	real_n = args.parallelize_n
nlines_per_job = str(math.ceil(sites_ct / real_n))

#split sites file into x files of nlines_per_job sites
subprocess.run(["split", sites_file, "-l", nlines_per_job, paralog_dir + "subsites"])

#prepare the array input script
subsites_files = glob.glob(paralog_dir + "subsites*")
paralog_array_input = scripts_dir + prefix + "_paralogARRAY_input.txt"
with open(paralog_array_input, 'w') as i:
	iterator = 1
	for ssf in subsites_files:
		i.write(str(iterator) + ":" + ssf + "\n")
		iterator += 1

#prepare the script to run mpileup + ngsParalog
paralog_script = scripts_dir + prefix + "_paralogARRAY.sh"
with open(paralog_script, 'w') as p:
	p.write("#!/bin/bash\n\n")
	p.write("#SBATCH --cpus-per-task=5\n")
	p.write("#SBATCH --time=1-00:00:00\n")
	p.write("#SBATCH --job-name=plog\n")
	p.write("#SBATCH --output=" + jobsout_dir + prefix + "_paralog_%A-%a.out\n")
	p.write("#SBATCH --error=" + jobsout_dir + prefix + "_paralog_%A-%a.err\n")
	p.write("#SBATCH --mail-type=FAIL\n")
	p.write("#SBATCH --mail-user=" + email + "\n")
	p.write("#SBATCH --array=1-" + str(int(args.parallelize_n)) + "%20\n\n")
	p.write("module unload bio/samtools/1.11\n")
	p.write("module load bio/samtools/1.11\n\n")

	p.write("JOBS_FILE=" + paralog_array_input + "\n")
	p.write("IDS=$(cat ${JOBS_FILE})\n\n")
	p.write("for sample_line in ${IDS}\n")
	p.write("do\n")
	p.write("""\tjob_index=$(echo ${sample_line} | awk -F ":" '{print $1}')\n""")
	p.write("""\tsites_file=$(echo ${sample_line} | awk -F ":" '{print $2}')\n""")
	p.write("\tif [[ ${SLURM_ARRAY_TASK_ID} == ${job_index} ]]; then\n")
	p.write("\t\tbreak\n")
	p.write("\tfi\n")
	p.write("done\n\n")

	p.write("samtools mpileup " + \
		"-b " + filtered_bamslist_file + " " + \
		"-l ${sites_file} " + \
		"-q 0 " + \
		"-Q 0 " + \
		"--ff UNMAP,DUP | " + \
		args.ngsParalog_executable + " calcLR " + \
		"-infile - " + \
		"-outfile ${sites_file}.lr " + \
		"-minQ 5 " + \
		"-minind 1 " + \
		"-mincov 1\n\n")

#RUN SIG TEST OVER CONCATENATED LR, INDEX SITES WITH ANGSD
wholegenome_lr_filename = paralog_dir + prefix + "_wholegenome.lr"
sigtest_lr_script = scripts_dir + prefix + "_sigLR.sh"
wgph_sites_filename = wholegenome_lr_filename.replace(".lr", "_retain.sites")
with open(sigtest_lr_script, 'w') as s:
	s.write("#!/bin/bash\n\n")
	s.write("#SBATCH --cpus-per-task=5\n")
	s.write("#SBATCH --time=0-01:00:00\n")
	s.write("#SBATCH --job-name=catLR\n")
	s.write("#SBATCH --output=" + jobsout_dir + prefix + "_concatLR_%A.out\n")
	s.write("#SBATCH --error=" + jobsout_dir + prefix + "_concatLR_%A.err\n")
	s.write("#SBATCH --mail-type=FAIL\n")
	s.write("#SBATCH --mail-user=" + email + "\n\n")

	s.write("module unload R/4.0.3 bio/angsd/0.933\n")
	s.write("module load R/4.0.3 bio/angsd/0.933\n\n")

	s.write("""echo -n "" > """ + wholegenome_lr_filename + "\n")
	s.write("cat " + '.lr '.join(subsites_files) + ".lr >> " + wholegenome_lr_filename + "\n\n")
	#use this instead of the cat line above if you find the chromchrom issue
#	s.write("for i in " + ' '.join(subsites_files) + "\n")
#	s.write("do\n")
#	s.write("LAST_LINE_LENGTH=$( tail -n 1 ${i} | tr -cd '\t' | wc -c | awk '{print $1+1}' )\n")
#	s.write("if (( ${LAST_LINE_LENGTH} == 5 )); then\n"
#	s.write("\tcat ${i} >> " + wholegenome_lr_filename + "\n")
#	s.write("else\n")
#	s.write("\tsed '$d' ${i} >> " + wholegenome_lr_filename + "\n")
#	s.write("fi\n")
#	s.write("done\n")
	s.write("Rscript --vanilla " + args.sigTest + " " + wholegenome_lr_filename + " " + args.alpha + "\n")
	s.write("angsd sites index " + wgph_sites_filename)

#GET GLS ACROSS EACH CHROM IN PARALLEL, TARGETING THE SITES THAT PASSED FILTERING
plmF_script = scripts_dir + prefix + "_polymorphic-filteredARRAY.sh"
with open(plmF_script, 'w') as plmF:
	plmF.write("#!/bin/bash\n\n")
	plmF.write("#SBATCH --cpus-per-task=10\n")
	plmF.write("#SBATCH --time=0-20:00:00\n")
	plmF.write("#SBATCH --job-name=plmF\n")
	plmF.write("#SBATCH --output=" + jobsout_dir + prefix + "_polymorphic-filtered_%A-%a.out\n")
	plmF.write("#SBATCH --error=" + jobsout_dir + prefix + "_polymorphic-filtered_%A-%a.err\n")
	plmF.write("#SBATCH --mail-type=FAIL\n")
	plmF.write("#SBATCH --mail-user=" + email + "\n")
	plmF.write("#SBATCH --array=1-" + str(len(chrs_list)) + "%24\n\n")

	plmF.write("module unload bio/angsd/0.933\n")
	plmF.write("module load bio/angsd/0.933\n\n")

	plmF.write("JOBS_FILE=" + scripts_dir + orig_prefix + "_angsdARRAY_input.txt\n")
	plmF.write("IDS=$(cat ${JOBS_FILE})\n\n")
	plmF.write("for sample_line in ${IDS}\n")
	plmF.write("do\n")
	plmF.write("""\tjob_index=$(echo ${sample_line} | awk -F ":" '{print $1}')\n""")
	plmF.write("""\tcontig=$(echo ${sample_line} | awk -F ":" '{print $2}')\n""")
	plmF.write("\tif [[ ${SLURM_ARRAY_TASK_ID} == ${job_index} ]]; then\n")
	plmF.write("\t\tbreak\n")
	plmF.write("\tfi\n")
	plmF.write("done\n\n")

	plmF.write("angsd -b " + filtered_bamslist_file + " " + \
		"-ref " + ref_genome + " " + \
		"-r ${contig}: " + \
		"-sites " + wgph_sites_filename + " " + \
		"-out " + gls_dir + prefix + "_${contig}_polymorphic-filtered " + \
		"-nThreads 10 " + \
		"-uniqueOnly 1 " + \
		"-remove_bads 1 " + \
		"-trim 0 " + \
		"-C 50 " + \
		"-doCounts 1 " + \
		"-doGlf 2 " + \
		"-GL 1 " + \
		"-doMaf 1 " + \
		"-doMajorMinor 1 " + \
		"-minMaf 0.05 " + \
		"-SNP_pval 1e-10 " + \
		"-minMapQ " + minq + " " + \
		"-minQ " + minq + " " + \
		"-setminDepth " + str(min_depth * filtered_n) + " " + \
		"-setmaxDepth " + str(max_depth * filtered_n) + " " + \
		"-doDepth 1 " + \
		"-dumpCounts 3")
	if endedness == "PE":
		plmF.write(" -only_proper_pairs 1")
	if min_ind is not None:
		plmF.write(" -minInd " + min_ind)

#LIST FILENAMES FOR CONCAT AND POSTERITY
plmF_beagle_files = []
plmF_maf_files = []
plmF_depths_files = []
plmF_counts_files = []

for chrom in chrs_list:
	plmF_beagle_file = gls_dir + prefix + "_" + chrom + "_polymorphic-filtered.beagle.gz"
	plmF_beagle_files.append(plmF_beagle_file)
	plmF_maf_file = gls_dir + prefix + "_" + chrom + "_polymorphic-filtered.mafs.gz"
	plmF_maf_files.append(plmF_maf_file)

#WHEN GLS CALC IS DONE, CONCAT BEAGLES AND MAFS, GET SITES FILE AND INDEX
#wgph beagle
all_chrs_wgph_beagle = gls_dir + prefix + "_wgph.beagle"
wgph_beagle_concat_script = scripts_dir + prefix + "_concatenate_wgph_beagles.sh"
with open(wgph_beagle_concat_script, 'w') as b:
	b.write("#!/bin/bash\n\n")
	b.write("#SBATCH --cpus-per-task=5\n")
	b.write("#SBATCH --job-name=cat-wgph-beagles\n")
	b.write("#SBATCH --mail-type=FAIL\n")
	b.write("#SBATCH --mail-user=" + email + "\n")
	b.write("#SBATCH --output=" + jobsout_dir + prefix + "_concatenate-wgph-beagles_%A.out\n")
	b.write("#SBATCH --error=" + jobsout_dir + prefix + "_concatenate-wgph-beagles_%A.err\n\n")

	b.write("zcat " + plmF_beagle_files[0] + " " + \
		"| head -n 1 > " + \
		all_chrs_wgph_beagle + "; " + \
		"for i in " + " ".join(plmF_beagle_files) + "\n" + \
		"do zcat $i " + \
		"| tail -n +2 -q >> " + \
		all_chrs_wgph_beagle + "\n" + \
		"done\n")
	b.write("gzip " + all_chrs_wgph_beagle)

#wgph mafs/sites
all_chrs_wgph_maf = gls_dir + prefix + "_wgph.mafs"
all_chrs_wgph_sites = gls_dir + prefix + "_wgph.sites"
wgph_mafs_concat_script = scripts_dir + prefix + "_concatenate_wgph_mafs.sh"
with open(wgph_mafs_concat_script, 'w') as m:
	m.write("#!/bin/bash\n\n")
	m.write("#SBATCH --cpus-per-task=5\n")
	m.write("#SBATCH --job-name=cat-wgph-mafs\n")
	m.write("#SBATCH --mail-type=FAIL\n")
	m.write("#SBATCH --mail-user=" + email + "\n")
	m.write("#SBATCH --output=" + jobsout_dir + prefix + "_concatenate-wgph-mafs_%A.out\n")
	m.write("#SBATCH --error=" + jobsout_dir + prefix + "_concatenate-wgph-mafs_%A.err\n\n")

	#in the event that a mafs file already exists, the following line will empty it (if the file doesn't exist, this line will initialize the file with nothing inside it)
	m.write("""echo -n "" > """ + all_chrs_wgph_maf + "\n")
	m.write("for i in " + " ".join(plmF_maf_files) + "\n" + \
		"do zcat $i | tail -n +2 -q >> " + \
		all_chrs_wgph_maf + "\n" + \
		"done\n")
	m.write("gzip " + all_chrs_wgph_maf + "\n\n")

	#if all went well, you can uncomment out the lines that remove subset files
#	c.write("for s in " + " ".join(subsites_files) + "\n")
#	c.write("do\n")
#	c.write("\trm $s\n")
#	c.write("done")

#Update ckpt file with everything
with open(args.ckpt_file, 'a') as ckpt:
	ckpt.write("paralogDIR\t" + paralog_dir + "\n")
	ckpt.write("glsFILES-polymorphic-filtered\t" + ",".join(plmF_beagle_files) + "\n")
	ckpt.write("mafsFILES-polymorphic-filtered\t" + ",".join(plmF_maf_files) + "\n")
	ckpt.write("wholegenomeBEAGLE-wgph\t" + all_chrs_wgph_beagle + ".gz\n")
	ckpt.write("wholegenomeMAF-wgph\t" + all_chrs_wgph_maf + ".gz\n")
	ckpt.write("wholegenomeSITES-wgph\t" + wgph_sites_filename + "\n")

#Print something helpful
print("Step 5 has finished successfully.")
print("Note that 'wgph' abbreviates 'whole genome polymorphic homologous', representing files of paralog-filtered (putatively homologous) SNPs across the whole genome.")
print("You will find a number of subset sites files, depending on the value you passed to parallelize_n.\n")
print("There are five new scripts:")
print("\t" + paralog_script + " generates an mpileup file for every subset sites file and runs ngsParalog.")
print("\t" + sigtest_lr_script + " concatenates the ngsParalog results from the subset sites files, identifies sites that are statistically likely to be homologous, and indexes those sites with angsd.")
print("\t" + plmF_script + " is very similar to step4, calculating genotype likelihoods for homologous sites, parallelized by chromosome.")
print("\t" + wgph_beagle_concat_script + " concatenates the beagle files resulting from executing " + plmF_script + ".")
print("\t" + wgph_mafs_concat_script + " concatenates the mafs files resulting from executing " + plmF_script + ".")
print("Both concatenation scripts can run simultaneously.\n")
print("Note that a whole-genome sites file is not generated from the whole-genome mafs file (as in step4), because the whole-genome sites file results from " + sigtest_lr_script + ".")

#import textwrap
#some_var = "the text to print out"
#print(textwrap.fill(some_var, [some value < 80]))
