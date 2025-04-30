#!/usr/bin/env python

import argparse
import os
import random
import subprocess
from collections import OrderedDict

#Read in config file for the run
parser = argparse.ArgumentParser()
parser.add_argument('--full_bamslist', '-f', help = 'The file listing all bamfiles (full paths included).')
parser.add_argument('--population_details', '-d', help = 'A comma-separated list of <population>:<number of individuals>.')
parser.add_argument('--population_pairs', '-p', help = 'A comma-separated list of <pop1>-<pop2>.')
parser.add_argument('--iteration', '-i', help = 'The iteration on which this script is being called.')
parser.add_argument('--sites_file', '-s', help = 'The .sites file of sites to include in the analysis.')
parser.add_argument('--reference_genome', '-r', help = 'The reference genome.')
parser.add_argument('--email', '-e', help = 'Email to notify about errors.')
parser.add_argument('--software_dir', '-a', help = 'Provide the path to the directory where your software is installed.')
parser.add_argument('--group_id', '-g', help = 'ID to differentiate parallel runs.')
args = parser.parse_args()

#make a directory to house all of this?
id_dir = args.group_id + "/"
if os.path.isdir(id_dir) is not True:
        os.mkdir(id_dir)

tmp_dir = id_dir + "tmp" + str(args.iteration) + "/"
if os.path.isdir(tmp_dir) is not True:
	os.mkdir(tmp_dir)

#parse info for each pop (name and sample size)
populations = OrderedDict()
pop_details_list = args.population_details.split(",")
for pop_and_ct in pop_details_list:
	populations[pop_and_ct.split(":")[0]] = [int(pop_and_ct.split(":")[1])]

#get set of bamfiles (I don't use sets very often, but it is important that every element is unique)
bams = set()
with open(args.full_bamslist, 'r') as b:
	for raw_bamfile_name in b:
		bams.add(raw_bamfile_name.rstrip())

#distribute bamfiles into pseudo-pops (sample size is equal to pop's true sample size)
for pop_id, n in populations.items():
	new_pop_bams = set(random.sample(list(bams), n[0]))
	pop_bamsfilename = tmp_dir + pop_id + "-" + args.iteration + "_bamslist.txt"
	n.append(pop_bamsfilename)
	out_filename = tmp_dir + pop_id + "-" + args.iteration
	n.append(out_filename)
	with open(pop_bamsfilename, 'w') as pb:
		for pop_bam in sorted(new_pop_bams):
			pb.write(pop_bam + "\n")
	bams -= new_pop_bams

#get the list of results files (each iteration will append to these files - a file for each population pair)
poppairs = OrderedDict()
for pop_pair in args.population_pairs.split(","):
	pop1 = pop_pair.split("-")[0]
	pop2 = pop_pair.split("-")[1]
	pop_pair_results_file = pop_pair + "_posterior.summary-fst.txt"
	poppairs[pop_pair] = [pop1, pop2, pop_pair_results_file]

#put all of this into a slurm script and sbatch it
scriptname = tmp_dir + "fst_iteration_" + args.iteration + ".sh"
with open (scriptname, 'w') as s:
	s.write("#!/bin/bash\n\n")
	s.write("#SBATCH --partition=bio\n")
	s.write("#SBATCH --nodes=1\n")
	s.write("#SBATCH --ntasks-per-node=24\n")
	s.write("#SBATCH --ntasks=24\n")
	s.write("#SBATCH --mem=214G\n")
	s.write("#SBATCH --job-name=" + str(args.iteration) + "\n")
	s.write("#SBATCH --mail-type=FAIL\n")
	s.write("#SBATCH --mail-user=" + args.email + "\n")
	s.write("#SBATCH --output=" + tmp_dir + "iteration" + str(args.iteration) + "_%A.out\n\n")

	s.write("ulimit -s unlimited\n")
	s.write("ulimit -l unlimited\n\n")

	s.write("module purge\n")
	s.write("module load bzip2/1.0.8 GCCcore/11.3.0\n\n")

	s.write("PATH=$PATH:" + args.software_dir + "angsd:" + args.software_dir + "angsd/misc\n\n")

	s.write("#Part 1#\n")
	for group, n_and_filename in populations.items():
		sample_ct = n_and_filename[0]
		bamsfilename = n_and_filename[1]
		outfilename = n_and_filename[2]
		s.write("###" + group + "###\n")
		s.write("angsd " + \
			"-b " + bamsfilename + " " + \
			"-ref " + args.reference_genome + " " + \
			"-anc " + args.reference_genome + " " + \
			"-sites " + args.sites_file + " " + \
			"-out " + outfilename + " " + \
			"-nThreads 20 " + \
			"-uniqueOnly 1 " + \
			"-remove_bads 1 " + \
			"-trim 0 " + \
			"-C 50 " + \
			"-minMapQ 15 " + \
			"-doCounts 1 " + \
			"-setminDepth " + str(sample_ct) + " " + \
			"-setmaxDepth " + str(float(sample_ct * 20)) + " " + \
			"-GL 1 " + \
			"-doGlf 1 " + \
			"-doMaf 1 " + \
			"-minMaf 0.05 " + \
			"-SNP_pval 1e-10 " + \
			"-doMajorMinor 1 " + \
			"-dumpCounts 3 " + \
			"-doDepth 1 " + \
			"-doSaf 1 " + \
			"-only_proper_pairs 1\n")
		s.write("realSFS " + \
			outfilename + ".saf.idx " + \
			"> " + outfilename + ".sfs\n\n")

	s.write("\n\n#Part 2#\n")
	for group_pair, accessories in poppairs.items():
		group1 = accessories[0]
		group1_saf = tmp_dir + group1 + "-" + args.iteration + ".saf.idx"
		group2 = accessories[1]
		group2_saf = tmp_dir + group2 + "-" + args.iteration + ".saf.idx"
		outfile_prefix = tmp_dir + group_pair + "-" + args.iteration
		results_file = id_dir + accessories[2]
		s.write("touch " + results_file + "\n")
		s.write("realSFS " + \
			group1_saf + " " + \
			group2_saf + " " + \
			"-P 20 " + \
			"-maxIter 30 " + \
			"> " + outfile_prefix + ".ml\n")
		s.write("realSFS fst index " + \
			group1_saf + " " + \
			group2_saf + " " + \
			"-sfs " + outfile_prefix + ".ml " + \
			"-fstout " + outfile_prefix + "\n")
		s.write("realSFS fst stats " + \
			outfile_prefix + ".fst.idx " + \
			"> " + outfile_prefix + ".summary.fst\n")
		s.write("head " + outfile_prefix + ".summary.fst >> " + results_file + "\n\n")

#	s.write("\n\n#Part 3#\n")
#	s.write("rm -r " + tmp_dir + "\n")

#subprocess.run(["sbatch", scriptname])

#if time is really an issue here
#step1_result = subprocess.run(["sbatch", array_script_step1])
#if step1_result.returncode == 0:
	#step2_result = subprocess.run(["sbatch", array_script_step2])
