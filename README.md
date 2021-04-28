# WGSfqs-to-genolikelihoods
A pipeline to assemble raw fastq files generated with whole-genome sequencing and calculate genotype likelihoods (on a SLURM scheduled HPC).

[built in python 3.6.8]
required packages (all standard): argparse, collections (OrderedDict), math, os, and subprocess

This pipeline exists to process raw lcWGS data and generates datasets based on 1) genotype-likelihoods, 2) called genotypes, and 3) SNPs. Currently, it only operates on SLURM-managed HPC clusters (but feel free to harvest the associated scripts in any way that facilitates your science). It was developed while consulting Dr. Nina Therkildsen's lcWGS Guide and Tutorial, available on github.

This guide is meant to be user-friendly as possible. Every time you pick it up, whether you're a first time user, you use it once a month, or you haven't looked at it in years, this manual should help you run with your data quickly. If something is confusing, please let Laura.Timm@noaa.gov know. If something included in the manual is insultingly over-explained, congratulations! You will likely be able to not only use this guide, but dig into the code and tailor it to your purposes.

Place the lcWGSpipeline scripts into your bin and make them executable before continuing.

Begin by making a directory that contains your gzipped fastq files. If you need to concatenate data, I recommend you try "concat_fqs.py", which is available at letimm's github. In addition to gunzipping, concatenating, and gzipping the concatenated file, this outputs a 'fastqs.txt' file, which lists all the gzipped fastq files in the directory. If you do not use "concat_fqs.py" to generate this list, you must provide it yourself. For the pipeline to name files consistently and correctly, file names MUST follow the format "sampleID_R1.fq.gz" (e.g., AGUL1_R1.fq.gz, ABLG645_R1.fq.gz, etc). For PE data, "R1" can change to "R2" (e.g., AGUL1_R1.fq.gz and AGUL1_R2.fq.gz, ABLG645_R1.fq.gz and ABLG645_R2.fq.gz, etc). The last thing needed in this directory is the config file. If you cloned this from letimm's github, this config file is named "lsWGS_config.txt" and it consists of a series of questions, separated from the user's responses by a TAB.
To recap, you want three things in this working directory {example file names, which will be used throughout this guide, are given in curly braces}:
1) the .fq.gz files {AGUL1_R1.fq.gz and ILIA52_R1.fq.gz}
2) the file listing the .fq.gz files {fastqs.txt}
3) and the config file that details how the pipeline will run {lcWGS_config.txt}.

This pipeline was developed on (largely) non-overlapping PE data, which means we treat it as SE data. This requires the PE files to be unzipped, concatenated, and deleted before gzipping the concatenated fq file. This can be accomplished with the python script "concat_fqs.py" (because this script globs every file ending with ".fq.gz" in the directory in which it is called, no input is needed). This script creates and submits a slurm script for each unique id (ABLG3_*_1.fq.gz and ABLG3_*_2.fq.gz have the uniqe id "ABLG3").

STEP0-CONFIGURE
The first step of the pipeline is to configure the SLURM space by loading modules, checking inputs, and indexing the reference genome (bwa-index).
Begin by running 'lcWGSpipeline_step0-configure.py -c <config_filename>'
This command takes one flag: '-c' or '--config_file'. This file must follow the format set out in lcWGS_config.txt; in fact, I recommend saving a copy of lcWGS_config.txt under some informative file name and editing it to reflect the pertinent information for your set-up.
After this step has run, you will have a new checkpoint file named according to the preferred prefix you provided (exampleprefix.ckpt).

STEP1-FASTQC
The next step takes the checkpoint file created in the last step and runs every fastq file through fastQC. Because many people complete this step on an HPC, this is accomplished by generating and executing a SLURM submission script for every fastq file. 
To run this step, call lcWGSpipeline_step1-fastqc.py -p <checkpoint_filename>'
This step does not alter the checkpoint file at all. 
It is recommended that you download the fastQC .html files to your local machine and open them in an internet brower to check for overall raw data quality before continuing. This can be accomplished with scp on Unix-based systems or pscp on MS-DOS/Windows machines. From your home machine, run 'scp <username>@<cluster_address>:/home/<username>/<working_dir>/fastqc/*.html Desktop'
This transfers fastQC files to the Desktop on your local machine.

STEP2-BWA
This step aligns reads to a reference genome, which was indexed in step 0. It adds the output .sam file names to the checkpoint file.
It is run with the command 'lcWGSpipeline_step2-bwa.py -p <checkpoint_filename>'

STEP3-SAMBAM
The SAM step runs the aligned reads through SAMTOOLS: 'fixmate' cleans up the read pairings and flags from BWA; a pair of 'view' statements converts the .sam file to a .bam file and filters the .bam file for non-unique and poor quality mappings; and 'sort' sorts the read pairings by coordinate (instead of read name). After a .bam file is built, duplicate reads are removed and (if the data is PE) overlapping reads are clipped to generate the final .bam. This step concludes by generating a new file (<prefix>_bams_list.txt), which lists all final .bam files that will serve as input to ANGSD. The checkpoint file is updated with the the name of the file containing the list of .bams. 
Run this step with 'lcWGSpipeline_step3-sambam.py -p <checkpoint_filename>'

STEP4-GLS
This is the heavy-lifting step, which filters the .bams according to user-specified flags and calculates genotype likelihoods before generating additional data types, specifically allele frequencies (a), genotypes (g), and snps (s). You may be interested in branching the pipeline at this step with the "-b" flag, specifying a file of branching parameterizations for filtering: 'lcWGSpipeline_step4-gls.py -p <checkpoint_filename> -b <branch_scheme_file> -t <data_output_types>'
The branch schema file must include at least one line: <branch_name>[TAB]<arguments to pass>. I suggest you use lcWGS-step4_branching_schema.txt, saved under an informative file name. The three parameterizations in lcWGS-step4_branching_schema.txt serve as useful starting places (although, if you have PE data, add "-only_proper_pairs 1" to each line) and lines can be added as needed. lcWGSpipeline_step4-filterdata generates and submits ANGSD commands to the HPC, so any line in the branch schema file will be added like so (in the example below, see arguments found within []):
If the branch schema file reads:
basic_filter	-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1
basic_filter_downsampled0.50	-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -c 50 -baq1 -downSample 0.50
lcWGSpipeline_step4-filterdata will submit the following to the cluster:
angsd -b <list_of_bamfiles_from_ckpt_file> -ref <reference_genome> -out <angsd_dir/[basic_filter]]> [-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1] -doCounts 1 -GL 2 -doGlf 4
angsd -b <list_of_bamfiles_from_ckpt_file> -ref <reference_genome> -out <angsd_dir/[basic_filter_downsampled0.50]> [-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -c 50 -baq1 -downSample 0.50] -doCounts 1 -GL 1 -doGlf 1
Some explanation of the flags added to the end (default bahviors)
+ genotype likelihoods (GLs) are calculated with the algorithm used in SAMtools [-GL 1] (NT recommends this when data uncertainty is high)
+ GLs are output as binary [-doGlf 1]
Should you choose to try new parameterizations after running one, remove any lines in the branch schema file that have already run (this prevents repetitive runs and wasting computational resources).
The checkpoint file is updated with a list of the genotype likelihood files (one for each branching scheme).
After these have been calculated, the default behaviors currently implemented in this output step are detailed below (in future, I may revise this step [the associated flag is given]):
ALLELE FREQUENCIES
+ major and minor alleles are fixed [-doMaf 1]
+ posterior probabilities are calculated using a uniform prior [-doPost 2] (instead of a prior informed by frequencies, 1)
+ both alleles are inferred from genotype likelihoods [-doMajorMinor 1]
GENOTYPES
+ genotypes (GTs) are called as major/minor alleles [-doGeno 9] (GTs printed as major/minor (1) + calculate posterior probabilities of all genotypes (8))
+ posterior probabilities use a uniform prior [-doPost 2]
+ genotypes will only be called if their posterior probability is greater than 0.75, otherwise it will be coded as missing data
SNPS
+ snps are only called if they have a minor allele frequency above 0.05 [-minMaf 0.05] and a p-value below 0.10 [-SNP-pval 0.10]
***I might want to enable another branching point here; maybe ask Wes if he has preferred limits***
