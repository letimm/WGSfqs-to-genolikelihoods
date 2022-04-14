# lcWGS Pipeline [built on python 3.6.8]
required packages (all standard): argparse, collections (OrderedDict), math, os, and subprocess

This pipeline exists to process raw lcWGS data and generates genotype likelihoods. Currently, it only operates on SLURM-managed HPC clusters (but feel free to harvest the associated scripts in any way that facilitates your science). It was developed while consulting Dr. Nina Therkildsen's lcWGS Guide and Tutorial, available on github.

This guide is meant to be as user-friendly as possible. Every time you pick it up, whether you're a first time user, you use it once a month, or you haven't looked at it in years, this manual should help you run with your data quickly. If something is confusing, please let Laura.Timm@noaa.gov know. If something included in the manual is insultingly over-explained, congratulations! You will likely be able to not only use this guide, but dig into the code and tailor it to your purposes.

Place the lcWGSpipeline scripts into your bin and make them executable before continuing. The majority of these scripts are in python, but a couple of R scripts also need to be in your bin and executable (HMM_log10FST+1_3norm.R [sourced from https://github.com/marqueda/HMM-detection-of-genomic-islands] and cov_to_eigen.R [slightly altered from source https://github.com/therkildsen-lab/genomic-data-analysis/blob/master/scripts/local_pca_2.R].

Begin by making a file that lists your gzipped fastq files (including the path). You must also prepare a config file. If you cloned this from letimm's github, this config file is named "lcWGS_config.txt" and it consists of a series of questions, separated from the user's responses by a TAB.
To recap, you want four things on the cluster {example file names, which will be used throughout this guide, are given in curly braces}:
1) the .fq.gz files {AGUL1_R1.fq.gz, AGUL1_R2.fq.gz and ILIA52_R1.fq.gz, ILIA52_R2.fq.gz} (note that data does not need to be paired-end),
2) the file listing the .fq.gz files {fastqs.txt},
3) the reference genome in fasta format (uncompressed),
4) and the config file that details how the pipeline will run {lcWGS_config.txt}.

# DATA ASSEMBLY
STEP0-CONFIGURE
command: lcWGSpipeline_step0-configure_vX.Y.py -c {lcWGS_config.txt}
flag(s): -c or --config_file
description: The first step of the pipeline is to configure the SLURM space by checking inputs, generating directories to house scripts and job outfiles, initializing a checkpoint file that will keep track of all the relevant information moving forward, and indexing the reference genome (bwa-index, samtools-fai).
output: After this step has run, you will have a new checkpoint file named according to the preferred prefix you provided {EXAMPLE.ckpt} as well as two scripts related to preparing the reference genome for analysis: {bwa-index_script.sh} and {fai_script.sh}.
success message: If this step runs successfully, you will see "Step 0 has finished successfully! You will find two new scripts in ./scripts/: {bwa-index_script.sh} and {fai_script.sh}. Both scripts can run simultaneously with 'sbatch'. However, there is no need to wait for these jobs to finish to move to step 1. To continue on, call step 1 to generate scripts for fastQC and multiQC. Remember to pass the newly-made checkpoint (.ckpt) file with the '-p' flag."
results: n/a

STEP1-QC
command: lcWGSpipeline_step1-fastqc_vX.Y.py -p {EXAMPLE.ckpt}
flag(s): -p or --ckpt_file
description: The next step takes the checkpoint file created in the last step and prepares scripts to quality-check the raw fastqs. I recommended that you download the multiQC .html file to your local machine and open it in an internet brower to check for overall raw data quality.
output: After this step has run, you will have two new scripts: {EXAMPLE-raw_fastqcARRAY.sh} (this file requires an array input file {EXAMPLE-raw_fqcARRAY_input.txt}, which is called by the fastqc script) and {EXAMPLE-raw_multiqcSLURM.sh}.
success message: "Step 1 has finished successfully! You will find two new scripts in ./scripts/: {EXAMPLE-raw_fastqcARRAY.sh} and {EXAMPLE-raw_multiqcSLURM.sh}. {EXAMPLE-raw_fastqcARRAY.sh} must run prior to launching {EXAMPLE-raw_multiqcSLURM.sh}. However, there is no need to wait for these jobs to finish to move to step 2. To continue on, call step 2 to generate a script for TRIMMOMATIC. Remember to pass the checkpoint (.ckpt) file with the '-p' flag."
results: multiqc.html

STEP2-TRIM
command: lcWGSpipeline_step2-trim_vX.Y.py -p {EXAMPLE.ckpt}
flag(s): -p or --ckpt_file
description: The next step takes the checkpoint file created in the last step and prepares scripts to trim adapters from raw fastqs and quality-check the trimmed fastqs. I recommended that you download the multiQC .html file to your local machine and open it in an internet brower to check for overall data quality before continuing to alignment.
output: After this step has run, you will have three new scripts: {EXAMPLE_trimARRAY.sh} (with array input file {EXAMPLE_trimARRAY_input.txt}), {EXAMPLE-trim_fastqcARRAY.sh} (with array input file {EXAMPLE-trim_fqcARRAY_input.txt}), and {EXAMPLE-trim_multiqcSLURM.sh}.
success message: "Step 2 has finished successfully! You will find three new scripts in ./scripts/: {EXAMPLE_trimARRAY.sh}, {EXAMPLE-trim_fastqcARRAY.sh}, {EXAMPLE-trim_multiqcSLURM.sh}. Each script must run in the above order and each job must finish before submitting the next. While there is no need to wait before running step 3 to generate the alignment script, it is probably wise to wait until the multiQC script has completed and the results have been viewed in a web browser prior to submitting the script written by step 3. Remember to pass the checkpoint (.ckpt) file with the '-p' flag."
results: {./trimmed/EXAMPLE.fq.gz}; multiqc1.html

STEP3-ALIGN
command: lcWGSpipeline_step3-align_vX.Y.py -p {EXAMPLE.ckpt}
flag(s): -p or --ckpt_file
description: This step aligns reads to a reference genome with BWA and runs the aligned reads through SAMTOOLS: 'fixmate' cleans up the read pairings and flags from BWA; a pair of 'view' statements converts the .sam file to a .bam file and filters the .bam file for non-unique and poor quality mappings; and 'sort' sorts the read pairings by coordinate (instead of read name). After a .bam file is built, duplicate reads are removed and (if the data is PE) overlapping reads are clipped to generate the final .bam. This step concludes by creating a file of depths for each alignment and generating a new file (<prefix>_bamslist.txt), which lists all final .bam files that will serve as input to ANGSD.
output: {EXAMPLE_alignARRAY.sh} (with array input file {EXAMPLE_alignARRAY_input.txt}) and {EXAMPLE_depthsARRAY.sh} (with array input file {EXAMPLE_depthsARRAY_input.txt}), as well as list of all bams ({EXAMPLE_bamslist.txt}).
success message: "Step 3 has finished successfully! You will find two new scripts in ./scripts/: {EXAMPLE_alignARRAY.sh} and {EXAMPLE_depthsARRAY.sh}. These scripts must be run in the order given above and both must finish before moving on to step 4. After {EXAMPLE_depthsARRAY.sh} has run, download {EXAMPLE_depths.csv} and generate a barchart of mean depths by individual. This will help determine whether any individuals fall substantially below the average depth (usually, we use a cutoff of 1x). If you identify samples with coverage that is 'too low', add the sample id(s) to a new file, referred to as the 'blacklist' of individuals to be excluded from genotype likelhiood calculations and the final data sets. After generating this blacklist, you can continue to step 4 to write scripts for generating the final data sets. Remember to pass the checkpoint (.ckpt) file with the '-p' flag AND the blacklist file with the '-b' flag."
results:

STEP4-DATA
command: lcWGSpipeline_step3-align_vX.Y.py -p {EXAMPLE.ckpt} -b {EXAMPLE_blacklist_1x.txt}
flag(s): -p or --ckpt_file; -b or --blacklist_inds
description: This step calculates genotype likelihoods for putatively polymorphic sites and all sites. For more detail and/or to tune parameters, open and revise the {EXAMPLE_globalARRAY.sh} and {EXAMPLE_polmorphicARRAY.sh} scripts.
output: {EXAMPLE_globalARRAY.sh} and {EXAMPLE_polymorphicARRAY.sh} (both of which use array input {EXAMPLE_angsdARRAY_input.txt}), as well as a list of the non-blacklisted bams ({filtered_bamslist.txt}).
success message: "Step 4 has finished successfully! You will find two new scripts in ./scripts/: {EXAMPLE_globalARRAY.sh} calculates genotype likelihoods across all sites on each chromosome (separately). {EXAMPLE_polymorphicARRAY.sh} calculates genotype likelihoods across all polymorphic sites on each chromosome (separately). Both scripts can run simultaneously. After they have run, you will have genotype likelihoods (gls) and allele frequencies (maf) for all sites in the genome (global) and putatively variable sites (polymorphic). Congratulations! This represents the end of data assembly."
results:

# DATA ANALYSIS
PCA & ADMIXTURE - revise me!
command: lcWGSpipeline_pca-admixture_vX.Y.py -p {EXAMPLE.ckpt} -k {10}
flag(s): -p or --ckpt_file; -k or --k_val_max (default = 10)
description: Running this step results in a script to generate a new beagle file (a concatenation of all genotype likelihoods across all chromosomes), which will serve as input for scripts to execute PCA and ADMIXTURE analysis. This step also generates a script to run PCA over each chromosome.
output: {EXAMPLE_concatenate_beagles.sh}, {EXAMPLE_pcangsdARRAY.sh} (which uses array input {EXAMPLE_angsdARRAY_input.txt}), and {EXAMPLE_wholegenome_pcangsd-admixture.sh}.
success message: "Three scripts have been generated: {EXAMPLE_concatenate_beagles.sh}, {EXAMPLE_pcangsdARRAY.sh}, and {EXAMPLE_wholegenome_pcangsd-admixture.sh}. {EXAMPLE_concatenate_beagles.sh} must run first. After it has finished, both {EXAMPLE_pcangsdARRAY.sh} and {EXAMPLE_wholegenome_pcangsd-admixture.sh} can run simultaneously. {EXAMPLE_pcangsdARRAY.sh} will run pca for each chromosome. {EXAMPLE_wholegenome_pcangsd-admixture.sh} will execute pca and admixture analysis for the whole genome polymorphic data."
results:

POPULATION-LEVEL ANALYSES
command: lcWGSpipeline_pop-analyses_vX.Y.py -p {EXAMPLE.ckpt} -g {pops.txt}
flag(s): -p or --ckpt_file; -g or --group_file (formatted as individual_ID<TAB>population)
description: This step generates a script to calculate a number of population-pair metrics. Site allele frequencies are calculated for each population, serving as input for the other calculations. Metrics include FST, dxy, and expected heterozygosity.
output: {EXAMPLE_popARRAY.sh} (which uses array input {EXAMPLE_angsdARRAY_input.txt}).
success message: "A new script has been generated to calculate population-level and populaiton-pair-level metrics: {EXAMPLE_popARRAY.sh}."
results:

HMM
command: lcWGSpipeline_hmm-bychrom_vX.Y.py -p {EXAMPLE.ckpt} -r {/home/ltimm/bin/HMM_log10FST+1_3norm.R}
flag(s): -p or --ckpt_file; -r or --r_script_loc (full path to the R script HMM_log10FST+1_3norm.R)
description: This step generates an array script which uses a Hidden Markov Model (HMM) to classify sites by fst value (high, low, background). This step runs by-chromosome over each population pair.
output: {EXAMPLE_hmmARRAY.sh} (which uses array input {EXAMPLE_hmmARRAY_input.txt}).
success message: "A new script has been generated to run HMM in R. Submitting {EXAMPLE_hmmARRAY.sh} will classify each site as high (1), low (2), or background differentiation (3), based on the fst value.
results:
Note: for this step to run correctly, you must have HMM_log10FST+1_3norm.R available for the script to access (in your path)."

ADAPTIVE REGIONS
command:
flag(s):
description:
output:
success message:
results:

LOSTRUCT
command:
flag(s):
description:
output:
success message:
results:
Note:

LINKAGE
command:
flag(s):
description:
output:
success message:
results:

GENOTYPE HEATMAPS
command:
flag(s):
description:
output:
success message:
results:
