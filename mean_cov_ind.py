#!/usr/bin/python3
import argparse
import gzip
import glob
import os

parser = argparse.ArgumentParser()
parser.add_argument('--infile', '-i', help = 'The name of the input file (a depth.gz file).')
parser.add_argument('--outfile', '-o', help = 'The name of the (comma-separated) outfile.')
args = parser.parse_args()

summed_depths = 0
num_depths = 0
mean_depth = 0
with gzip.open(args.infile, 'rt') as depths_in:
		for line in depths_in:
				depth = float(line.rstrip())
				summed_depths += depth
				num_depths += 1
		mean_depth = summed_depths / num_depths

#if os.path.isfile(args.outfile):
with open(args.outfile, 'a') as d:
	d.write(args.infile + "\t" + str(mean_depth) + "\n")
#else:
#	with open(args.outfile, 'w') as w:
#		w.write(args.infile + "\t" + str(mean_depth) + "\n")
