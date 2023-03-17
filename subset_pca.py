#!/usr/bin/python3
import argparse
import gzip

parser = argparse.ArgumentParser()
parser.add_argument('--beagle_file', '-i', help = 'a .beagle.gz file (the "global.beagle.gz" file from ANGSD)')
parser.add_argument('--fsts_file', '-f', required = False, help = 'a file of fst files, with each fst file on a new line (only required if a max and/or min fst are set)')
parser.add_argument('--chromosome', '-c', help = 'the chromosome on which the region/island/peak occurs (as labeled in the reference genome fasta')
parser.add_argument('--prefix', '-p', help = 'a unique identifier to differentiate the run')
parser.add_argument('--start_pos', '-s', type = int, help = 'start position of the region of interest')
parser.add_argument('--end_pos', '-e', type = int, help = 'end position of the region of interest')
parser.add_argument('--max_fst', '-m', type = float, help = 'maximum fst of interest', default = 1.0)
parser.add_argument('--min_fst', '-l', type = float, help = 'minimum fst of interest', default = 0.0)
parser.add_argument('--bamslist', '-b', help = 'file including a list of bamfiles to include in pca')
parser.add_argument('--ref_genome', '-r', help = 'reference genome')
args = parser.parse_args()

sites_list = []
if args.max_fst != 1.0 or args.min_fst != 0.0:
	fsts_filelist = []
	with open(args.fsts_file, 'r') as f:
		for filename in f:
			fsts_filelist.append(filename.rstrip())
	for fst_file in fsts_filelist:
		with open(fst_file, 'r') as fst:
			next(fst)
			for fst_line in fst:
				chrom = fst_line.rstrip().split('\t')[1]
				if chrom == args.chromosome:
					site = int(fst_line.rstrip().split('\t')[2])
					if args.start_pos <= site <= args.end_pos:
						fst_val = float(fst_line.rstrip().split('\t')[4])
						if fst_val < 0:
							fst_val = 0
						if args.min_fst <= fst_val <= args.max_fst:
							if site not in sites_list:
								sites_list.append(site)
	sites_list.sort()
else:
	i = args.start_pos
	while i <= args.end_pos:
		sites_list.append(i)
		i += 1

subset_beagle_outfile = args.prefix + ".beagle.gz"
with gzip.open(args.beagle_file, 'rt') as b:
	linenum = 0
	with gzip.open(subset_beagle_outfile, 'wb') as newb:
		for site_line in b:
			if linenum == 0:
				newb.write(site_line.encode())
				linenum += 1
			else:
				pos = int(site_line.rstrip().split('\t')[0].split("_")[1])
				if min(sites_list) <= pos <= max(sites_list):
					if pos in sites_list:
						newb.write(site_line.encode())
				elif pos > max(sites_list):
					break

pca_script = args.prefix + ("_pcaSLURM.sh")
with open(pca_script, 'w') as s:
	s.write("#!/bin/bash\n\n")
	s.write("#SBATCH --cpus-per-task=5\n")
	s.write("#SBATCH --job-name=" + args.prefix + "_pca\n")
	s.write("#SBATCH --output=" + args.prefix + "_pca.out\n\n")
	s.write("module unload bio/pcangsd/0.99\n")
	s.write("module load bio/pcangsd/0.99\n")
	s.write("source /opt/bioinformatics/venv/pcangsd-0.99/bin/activate\n\n")
	s.write("pcangsd.py " + \
		"-threads 5 " + \
		"-beagle " + subset_beagle_outfile + " " + \
		"-o " + args.prefix + "_pca " + \
		"-sites_save " + \
		"-pcadapt")