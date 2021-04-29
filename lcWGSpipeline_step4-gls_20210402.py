#!/usr/bin/python3
import argparse
import math
import os
import subprocess
from collections import OrderedDict

#Read in config file for the run
parser = argparse.ArgumentParser()
parser.add_argument('--ckpt_file', '-p', help = 'Please provide the checkpoint file created in step 0.')
parser.add_argument('--branch_schema', '-b', help = 'Please provide a file containing parameterizations (branches) for the data filtering step.')
args = parser.parse_args()

#Initialize run config variables with some default values
scripts_dir = None
bwa_dir = None
bamtools_dir = None
ref_genome = None
angsd_dir = None
max_nodes = None
prefix = None
bams_listfile = None
branches = OrderedDict()
gls_filenames = []

#Parse the config file to determine what's needed (if user wants bwa, no need for picard or bowtie2)
with open(args.ckpt_file, 'r') as last_step_ckpt:
	for raw_ckpt_line in last_step_ckpt:
		ckpt_line = raw_ckpt_line.rstrip()
		ckpt_setting = ckpt_line.split('\t')
		if ckpt_setting[0] == "scriptsDIR":
			scripts_dir = ckpt_setting[1]
		elif ckpt_setting[0] == "bamtoolsDIR":
			bamtools_dir = ckpt_setting[1]
		elif ckpt_setting[0] == "refgenomeFASTA":
			ref_genome = ckpt_setting[1]
		elif ckpt_setting[0] == "angsdDIR":
			angsd_dir = ckpt_setting[1]
		elif ckpt_setting[0] == "maxNODES":
			max_nodes = float(ckpt_setting[1])
		elif ckpt_setting[0] == "prefix":
			prefix = ckpt_setting[1]
		elif ckpt_setting[0] == "bamsLIST":
			bams_listfile = ckpt_setting[1]

#Parse the parameterizations in the branching schema file
n_branches = 0
with open(args.branch_schema, 'r') as branchfile:
	for scheme in branchfile:
		raw_scheme_line = scheme.rstrip()
		label_scheme = raw_scheme_line.split('\t')
		branches[label_scheme[0]] = label_scheme[1]
		n_branches += 1

#Determine how many nodes can be used to run each job
n_nodes = math.floor(max_nodes / float(n_branches))
if max_nodes == 1:
	n_nodes = 1
elif n_nodes <= 2:
	n_nodes = 2
else:
	n_nodes = n_nodes

gls_dir = angsd_dir + "genotype_likelihoods/"
if os.path.isdir(gls_dir) is not True:
	os.mkdir(gls_dir)

for label, scheme in branches.items():
	gls_filename = prefix + "_" + label + ".glf.gz"
	gls_filenames.append(gls_filename)
	script = scripts_dir + prefix + "_" + label + "_glsSLURM.sh"
	with open(script, 'w') as s:
		s.write("#!/bin/bash\n\n")
		s.write("#SBATCH --nodes=" + str(n_nodes) + "\n")
		s.write("#SBATCH --ntasks=" + str(n_nodes) + "\n")
		s.write("#SBATCH --job-name=gls_" + label + "\n")
		s.write("#SBATCH --output=" + gls_dir + prefix + "_" + label + "_gls.out\n\n")
		s.write("module unload bio/angsd\n")
		s.write("module load bio/angsd\n\n")
		s.write("angsd -b " + bams_listfile + \
			" -ref " + ref_genome + \
			" -out " + gls_dir + prefix + "_" + label + " " + scheme + \
			" -GL 1 -doGlf 1 -doMajorMinor 1 -doMaf 1 -doGeno 1 -doPost 2\n")
	subprocess.call(["sbatch", script])

#Update checkpoint file with the genotype likelihood file names (one for each branching scheme)
with open(args.ckpt_file, 'a') as ckpt:
	ckpt.write("projNODES\t" + str(n_nodes) + "\n")
	for glsfile in gls_filenames:
		ckpt.write("GLS\t" + glsfile + "\n")