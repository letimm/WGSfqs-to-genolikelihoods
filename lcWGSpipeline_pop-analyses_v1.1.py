#!/usr/bin/python3
import argparse
import math
import os
from collections import OrderedDict

#Read in config file for the run
parser = argparse.ArgumentParser()
parser.add_argument('--ckpt_file', '-p', help = 'Please provide the checkpoint file created in step 0.')
parser.add_argument('--group_file', '-g', help = 'Please provide a tab-delimited file assigning each individual to a group (ind_id<TAB>group_assignment).')
args = parser.parse_args()

#Initialize run config variables with some default values
working_dir = None
bam = None
scripts_dir = None
jobsout_dir = None
gls_dir = None
prefix = None
email = None
ref_fasta = None
endedness = None
pop_bamfiles = OrderedDict()
chrs = []
polymorphic_sites_file = None

#Parse the ckpt file
with open(args.ckpt_file, 'r') as last_step_ckpt:
	for raw_ckpt_line in last_step_ckpt:
		ckpt_line = raw_ckpt_line.rstrip()
		ckpt_setting = ckpt_line.split('\t')
		if ckpt_setting[0] == "workingDIR":
			working_dir = ckpt_setting[1]
		elif ckpt_setting[0] == "bamsLIST-filtered":
			bam = ckpt_setting[1]
		elif ckpt_setting[0] == "scriptsDIR":
			scripts_dir = ckpt_setting[1]
		elif ckpt_setting[0] == "jobsoutDIR":
			jobsout_dir = ckpt_setting[1]
		elif ckpt_setting[0] == "glsDIR":
			gls_dir = ckpt_setting[1]
		elif ckpt_setting[0] == "prefix":
			prefix = ckpt_setting[1]
		elif ckpt_setting[0] == "email":
			email = ckpt_setting[1]
		elif ckpt_setting[0] == "refgenomeFASTA":
			ref_fasta = ckpt_setting[1]
		elif ckpt_setting[0] == "ENDEDNESS":
			endedness = ckpt_setting[1]
		elif ckpt_setting[0] == "chrsLIST":
			chrs = ckpt_setting[1].split(",")
		elif ckpt_setting[0] == "groupsBAMS":
			pop_list = ckpt_setting[1].split(",")
			for file_and_count in pop_list:
				pop_name = file_and_count.split(":")[0]
				bamfile = file_and_count.split(":")[1]
				count = file_and_count.split(":")[2]
				pop_bamfiles[pop_name] = [bamfile,count]
		elif ckpt_setting[0] == "wholegenomeSITES-polymorphic":
			polymorphic_sites_file = ckpt_setting[1]


div_dir = working_dir + "diversity/"
fst_dir = working_dir + "fst/"
if os.path.isdir(div_dir) is not True:
	os.mkdir(div_dir)
	with open(args.ckpt_file, 'a') as ckpt:
		ckpt.write("diversityDIR\t" + div_dir + "\n")
if os.path.isdir(fst_dir) is not True:
	os.mkdir(fst_dir)
	with open(args.ckpt_file, 'a') as ckpt:
		ckpt.write("fstDIR\t" + fst_dir + "\n")


if len(list(pop_bamfiles.keys())) == 0:
	bam_files = []
	with open(bam, 'r') as b:
		for raw_bam_line in b:
			bam_filename = raw_bam_line.rstrip()
			bam_files.append(bam_filename)

	groups = OrderedDict()
	with open(args.group_file, 'r') as grps:
		for raw_line in grps:
			ind = raw_line.rstrip().split('\t')[0] + "_"
			grp = raw_line.rstrip().split('\t')[1]
			for bf in bam_files:
				if ind in bf:
					ind_bamfile = bf
					if grp in groups.keys():
						groups[grp].append(bf)
					else:
						groups[grp] = [bf]

	pop_info_list = []
	for group, inds in groups.items():
		group_filename = working_dir + prefix + "_" + group + "_bams.txt"
		pop_bamfiles[group] = [group_filename]
		sample_ct = 0
		with open(group_filename, 'w') as gr_out:
			for sample_bam in inds:
				gr_out.write(sample_bam + '\n')
				sample_ct += 1
		pop_bamfiles[group].append(sample_ct)
		pop_info = group + ":" + group_filename + ":" + str(sample_ct)
		pop_info_list.append(pop_info)

	with open(args.ckpt_file, 'a') as ckpt:
		ckpt.write("groupsBAMS\t" + ",".join(pop_info_list) + "\n")

pops_list = list(pop_bamfiles.keys())

pop_script = scripts_dir + prefix + "_popARRAY.sh"
with open(pop_script, 'w') as pop:
	pop.write("#!/bin/bash\n\n")
	pop.write("#SBATCH --cpus-per-task=10\n")
	pop.write("#SBATCH --time=1-00:00:00\n")
	pop.write("#SBATCH --job-name=pop_" + prefix + "\n")
	pop.write("#SBATCH --output=" + jobsout_dir + prefix + "_pop-analyses_%A-%a.out\n")
	pop.write("#SBATCH --mail-type=FAIL\n")
	pop.write("#SBATCH --mail-user=" + email + "\n")
	pop.write("#SBATCH --array=1-" + str(len(chrs)) + "%24\n\n")
	pop.write("module unload bio/angsd/0.933 bio/ngstools/202202\n")
	pop.write("module load bio/angsd/0.933 bio/ngstools/202202\n\n")

	pop.write("JOBS_FILE=" + scripts_dir + prefix + "_angsdARRAY_input.txt\n")
	pop.write("IDS=$(cat ${JOBS_FILE})\n\n")
	pop.write("for sample_line in ${IDS}\n")
	pop.write("do\n")
	pop.write("""\tjob_index=$(echo ${sample_line} | awk -F ":" '{print $1}')\n""")
	pop.write("""\tchrom=$(echo ${sample_line} | awk -F ":" '{print $2}')\n""")
	pop.write("\tif [[ ${SLURM_ARRAY_TASK_ID} == ${job_index} ]]; then\n")
	pop.write("\t\tbreak\n")
	pop.write("\tfi\n")
	pop.write("done\n\n")

	for population in sorted(pops_list):
		pop.write("angsd -b " + pop_bamfiles[population][0] + " " + \
			"-r ${chrom}: " + \
			"-sites " + polymorphic_sites_file + " " + \
			"-ref " + ref_fasta + " " + \
			"-anc " + ref_fasta + " " + \
			"-out " + div_dir + prefix + "_${chrom}_" + population + "_polymorphic_folded " + \
			"-nThreads 10 " + \
			"-uniqueOnly 1 " + \
			"-remove_bads 1 " + \
			"-trim 0 " + \
			"-C 50 " + \
			"-minMapQ 15 " + \
			"-doCounts 1 " + \
			"-setminDepth " + str(pop_bamfiles[population][1]) + " " + \
			"-setmaxDepth " + str(float(pop_bamfiles[population][1]) * 20) + " " + \
			"-GL 1 " + \
			"-doGlf 3 " + \
			"-doMaf 1 " + \
			"-doMajorMinor 1 " + \
			"-doSaf 1")
		if endedness == "PE":
			pop.write(" -only_proper_pairs 1")
		pop.write("\n\n")

		#calculate site frequency spectra to get diversity and selection values (theta, Watterson, Tajima's D)
		#pop.write("### Calculate site frequency spectra to get diversity and selection values (theta, Watterson, Tajima's D) for each population. ###\n")
		#pop.write("realSFS " + div_dir + prefix + "_${chrom}_" + population + "_polymorphic_folded.saf.idx " + \
		#	"-fold 1 " + \
		#	"-P 10 > " + \
		#	div_dir + prefix + "_${chrom}_" + population + "_polymorphic_folded.saf.sfs\n")
		#pop.write("realSFS saf2theta " + div_dir + prefix + "_${chrom}_" + population + "_polymorphic_folded.saf.idx " + \
		#	"-fold 1 " + \
		#	"-sfs " + div_dir + prefix + "_${chrom}_" + population + "_polymorphic_folded.saf.sfs " + \
		#	"-outname " + div_dir + prefix + "_${chrom}_" + population + "_polymorphic_folded\n")
		#pop.write("thetaStat do_stat " + \
		#	div_dir + prefix + "_${chrom}_" + population + "_polymorphic_folded.thetas.idx " +\
		#	"-win 1 " + \
		#	"-step 1 " + \
		#	"-outnames " + div_dir + prefix + "_${chrom}_" + population + "_polymorphic_folded.thetasWindow\n")
		#pop.write("thetaStat print " + \
		#	div_dir + prefix + "_${chrom}_" + population + "_polymorphic_folded.thetas.idx > " + \
		#	div_dir	+ prefix + "_${chrom}_" + population + "_polymorphic_folded.thetas.txt\n\n")

		#inbreeding
		#pop.write("### Calculate inbreeding for each population. ###\n")
		#pop.write("SITES_CT=`zcat " + \
		#	div_dir	+ prefix + "_${chrom}_" + population + "_polymorphic_folded.mafs.gz " + \
		#	"| tail -n+2 | wc -l`\n")
		#pop.write("zcat " + \
		#	div_dir	+ prefix + "_${chrom}_" + population + "_polymorphic_folded.glf.gz > " + \
		#	div_dir + prefix + "_${chrom}_" + population + "_polymorphic_folded.glf\n")
		#pop.write("ngsF " + \
		#	"--n_ind " + str(pop_bamfiles[population][1]) + " " + \
		#	"--n_sites ${SITES_CT} " + \
		#	"--glf " + div_dir	+ prefix + "_${chrom}_" + population + "_polymorphic_folded.glf " + \
		#	"--out " + inbre_dir	+ prefix + "_${chrom}_" + population + "_polymorphic_folded.indF\n\n")

	#look at population pairs by chromosome
	processed_pops = []
	for pop1 in sorted(pops_list):
		processed_pops.append(pop1)
		for pop2 in sorted(pops_list):
			if pop2 not in processed_pops:
				pop.write("realSFS " + \
					div_dir + prefix + "_${chrom}_" + pop1 + "_polymorphic_folded.saf.idx -fold 1 " + \
					div_dir + prefix + "_${chrom}_" + pop2 + "_polymorphic_folded.saf.idx -fold 1 > " + \
					fst_dir + prefix + "_${chrom}_" + pop1 + "-" + pop2 + "_polymorphic_folded.sfs\n")
				pop.write("realSFS fst index " + \
					div_dir + prefix + "_${chrom}_" + pop1 + "_polymorphic_folded.saf.idx -fold 1 " + \
					div_dir + prefix + "_${chrom}_" + pop2 + "_polymorphic_folded.saf.idx -fold 1 " + \
					"-sfs " + fst_dir + prefix + "_${chrom}_" + pop1 + "-" + pop2 + "_polymorphic_folded.sfs " + \
					"-fstout " + fst_dir + prefix + "_${chrom}_" + pop1 + "-" + pop2 + "_polymorphic_folded.sfs.pbs " + \
					"-whichFst 1\n")
				pop.write("realSFS fst stats2 " + \
					fst_dir + prefix + "_${chrom}_" + pop1 + "-" + pop2 + "_polymorphic_folded.sfs.pbs.fst.idx " + \
					"-win 1 " + \
					"-step 1 > " + \
					fst_dir + prefix + "_${chrom}_" + pop1 + "-" + pop2 + "_polymorphic_folded.sfs.pbs.fst.txt\n\n")

				#identify intersecting sites for which to calculate summary stats (seg sites, He, dxy, fixed sites)
				#pop.write("### Calculate summary statistics for every population and population pair (segregating sites, expected heterozygosity, dxy, and fixed sites). ###\n")
				#pop.write("realSFS print " + \
				#	div_dir + prefix + "_${chrom}_" + pop1 + "_polymorphic_folded.saf.idx " + \
				#	div_dir + prefix + "_${chrom}_" + pop2 + "_polymorphic_folded.saf.idx " + \
				#	"| cut -f 1-2 > " + \
				#	dxy_dir + prefix + "_${chrom}_" + pop1 + "-" + pop2 + "_polymorphic_folded_intersecting-sites.txt\n")
				#pop.write("NSITES=`wc -l " + \
				#	dxy_dir + prefix + "_${chrom}_" + pop1 + "-" + pop2 + "_polymorphic_folded_intersecting-sites.txt " + \
				#	'| cut -f 1 -d " "`\n')
				#pop.write("zcat " + \
				#	div_dir + prefix + "_${chrom}_" + pop1 + "_polymorphic_folded.saf.gz > " + \
				#	dxy_dir + prefix + "_${chrom}_" + pop1 + "_polymorphic_folded.saf\n")
				#pop.write("zcat " + \
				#	div_dir + prefix + "_${chrom}_" + pop2 + "_polymorphic_folded.saf.gz > " + \
				#	dxy_dir + prefix + "_${chrom}_" + pop2 + "_polymorphic_folded.saf\n\n")
				#pop.write("ngsStat " + \
				#	"-npop 2 " + \
				#	"-postfiles " + \
				#	dxy_dir + prefix + "_${chrom}_" + pop1 + "_polymorphic_folded.saf " + \
				#	dxy_dir + prefix + "_${chrom}_" + pop2 + "_polymorphic_folded.saf " + \
				#	"-nsites ${NSITES} " + \
				#	"-nind " + \
				#	str(pop_bamfiles[pop1][1]) + " " + \
				#	str(pop_bamfiles[pop2][1]) + " " + \
				#	"-outfile " + \
				#	dxy_dir + prefix + "_${chrom}_" + pop1 + "-" + pop2 + "_polymorphic_folded_stats.txt\n")
				#pop.write("\n")

print("The resulting script calculates population-level and population-pair-level metrics:")
print(pop_script + " performs a genome scan and calculates fst at every polymorphic site across the genome.")
print("This is parallelized by chromosome.\n")
