#!/usr/bin/python3
import argparse
from scipy import stats
from collections import OrderedDict
from statistics import mean
from statistics import stdev

#Read in config file for the run
parser = argparse.ArgumentParser()
parser.add_argument('--posterior_distributions', '-d', help = 'The file listing all files of posteriors (*_posterior.summary-fst.txts output from generate_fst_posterior.py)')
parser.add_argument('--fsts', '-f', help = 'The file listing all files of fsts (*.global.fsts output from summary-fst_pt2).')
parser.add_argument('--prefix', '-p', help = 'The prefix for the outfile.', default = "fst_significance_test")
group = parser.add_mutually_exclusive_group()
group.add_argument('--weighted', action='store_true')
group.add_argument('--unweighted', action='store_true')
args = parser.parse_args()

#NOTE! You need to setup a virtual environment in which to execute this script.
#See https://docs.google.com/document/d/1nn0T0OWEsQCBoCdaH6DSY69lQSbK3XnPlseyyQuU2Lc/edit for directions, but if you want to just copy-paste...
#   python3 -m venv ~/bin/sciencesnake
#   source ~/bin/sciencesnake/bin/activate
#   pip install scipy
#   fst_significance_test.py -h
#when you're done, you can exit the venv with
#   deactivate

pair_fst = OrderedDict() #[pair_label] = fst_val
pair_dist = OrderedDict() #[pair_label] = [fst_distribution]

#set whether to record weighted or unweighted fst
result_index = None
if args.unweighted:
	result_index = 0
elif args.weighted:
	result_index = 1

#parse fst files
fst_files = []
with open(args.fsts, 'r') as fst_files_file:
	for fst_file in fst_files_file:
		fst_files.append(fst_file.rstrip())
for fst_f in fst_files:
	#parse filename for pop_pair label
	pair_label = fst_f.split("/")[-1].split("_")[-1].split(".")[0]
	#record data from file
	with open(fst_f, 'r') as f:
		for weighted_and_unweighted_fst_data_line in f:
			raw_fst = float(weighted_and_unweighted_fst_data_line.rstrip().split('\t')[result_index])
			real_fst = None
			if raw_fst < 0:
				real_fst = 0
			else:
				real_fst = raw_fst
			pair_fst[pair_label] = real_fst

#parse distribution files
dist_files = []
with open(args.posterior_distributions, 'r') as pd_files_file:
	for pd_file in pd_files_file:
		dist_files.append(pd_file.rstrip())
for pd_f in dist_files:
	#parse filename for pop_pair label
	pair_label = pd_f.split("/")[-1].split("_")[0]
	#record data from file
	dist_vals = []
	with open(pd_f, 'r') as p:
		for weighted_and_unweighted_fst_dist_line in p:
			raw_fst = float(weighted_and_unweighted_fst_dist_line.rstrip().split('\t')[result_index])
			if raw_fst < 0:
				real_fst = 0
			else:
				real_fst = raw_fst
			dist_vals.append(real_fst)
	pair_dist[pair_label] = dist_vals

#make a list of pop-pair keys (should be universal across dicts...)
pop_pairs = list(pair_fst.keys())

#use the list of keys to traverse the dicts and write an outfile
sig_results = OrderedDict()
for pair in sorted(pop_pairs):
	fst_result = pair_fst[pair]
	dist_mean = mean(pair_dist[pair])

#fsts seem to empirically follow an exponential distribution (see Elhaik 2012 [https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0049837])
	dist_lambda = 1 / dist_mean
	cdf_result = stats.expon(loc = dist_mean, scale = dist_mean).cdf(fst_result)
	pval = 1 - cdf_result
	sig_results[pair] = ["%.3f" % fst_result, "%.3f" % dist_mean, "%.3f" % dist_lambda, "%.3f" % cdf_result, "%.3f" % pval]

outfilename = args.prefix + "_sigtest.tsv"
with open(outfilename, 'w') as o:
	fst_header = None
	if args.weighted:
		fst_header = "weighted_fst_value"
	elif args.unweighted:
		fst_header = "unweighted_fst_value"
	o.write("group_pair\t" + fst_header + "\tdistribution_avg\tdistribution_lambda\texponential_cdf\tp-val\n")
	for ppair, results in sig_results.items():
		o.write(ppair + "\t" + "\t".join(results) + "\n")
