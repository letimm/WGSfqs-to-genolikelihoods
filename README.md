# RTFM: lcWGS Pipeline [built on python 3.6.8]

This pipeline exists to process raw lcWGS data and generates genotype likelihoods. Currently, it only operates on SLURM-managed HPC clusters (but feel free to harvest the associated scripts in any way that facilitates your science). It was developed while consulting Dr. Nina Therkildsen's lcWGS Guide and Tutorial, available on github.

This guide is meant to be as user-friendly as possible. Every time you pick it up, whether you're a first time user, you use it once a month, or you haven't looked at it in years, this manual should help you run with your data quickly. If something is confusing, please let Laura.Timm@noaa.gov know. If something included in the manual is insultingly over-explained, congratulations! You will likely be able to not only use this guide, but dig into the code and tailor it to your purposes.

If you feel that this pipeline contributes meaningfully to your project, please cite the version publicly available on LE Timm's github: https://github.com/letimm/WGSfqs-to-genolikelihoods

Place the lcWGSpipeline scripts into your bin and make them executable before continuing. The majority of these scripts are in python, but a couple of R scripts also need to be in your bin and executable (HMM_log10FST+1_3norm.R [sourced from https://github.com/marqueda/HMM-detection-of-genomic-islands] and cov_to_eigen.R [slightly altered from source https://github.com/therkildsen-lab/genomic-data-analysis/blob/master/scripts/local_pca_2.R].

Begin by making a file that lists your gzipped fastq files (including the path). You must also prepare a config file. If you cloned this from letimm's github, this config file is named "lcWGS_config.txt" and it consists of a series of questions, separated from the user's responses by a TAB. Some users have reported that TABs can become SPACEs during the transfer, so it's worth confirming whitespaces after "?s" are TABs.

To begin, you want a few things on the cluster (*example file names* are italicized):
1) the .fq.gz files (*ABLG100_R1.fq.gz*, *ABLG100_R2.fq.gz* and *ABLG101_R1.fq.gz*, *ABLG101_R2.fq.gz*; note that data does not need to be paired-end),
2) the file listing the .fq.gz files (*./example_files/fastqs.txt*),
3) the reference genome in fasta format (uncompressed) (*ref.fa*),
4) the file of adapters used to sequence the reads (*./example_files/NexteraPE-PE.fa*}
5) a file listing the chromosomes/contigs to analyze (*./example_files/chromosomes.txt*)
6) and the config file that details how the pipeline will run (*./example_files/lcWGS_config.txt*).

If you don't have immediate access to this README, details of how to run any python script can be accessed with
`scriptname.py -h`

### DATA ASSEMBLY  
##### STEP0 - CONFIGURE  
**command:** 
```
lcWGSpipeline_step0-configure.py -c lcWGS_config.txt
```  
**flag(s):** -c or --config_file  

**description:** The first step of the pipeline is to configure the SLURM space by checking inputs, generating directories to house scripts and job outfiles, initializing a checkpoint file that will keep track of all the relevant information moving forward, and indexing the reference genome (bwa-index, samtools-fai).
The checkpoint file that is created in this step is updated throughout the pipeline with relevant parameters at each step.\

**output:** 
1) After this step has run, you will have a new checkpoint file named according to the preferred prefix you provided `EXAMPLE.ckpt` as well as two scripts related to preparing the reference genome for analysis: 
2) `bwa-index_script.sh` and 
3) `fai_script.sh`.

**success message:** Step 0 has finished successfully! You will find two new scripts in `./scripts`: `ref_bwa-index.sh` and `ref_fai.sh`. Both scripts can run simultaneously with 'sbatch'. However, there is no need to wait for these jobs to finish to move to step 1. To continue on, call step 1 to generate scripts for fastQC and multiQC. Remember to pass the newly-made checkpoint (.ckpt) file with the '-p' flag.

**results:**
1) ref.fai (after running `ref_fai.sh`) 
after running `ref_bwa-index.sh`: 
2) ./bwa/ref.amb
3) ./bwa/ref.ann
4) ./bwa/ref.bwt
5) ./bwa/ref.pac
6) ./bwa/ref.sa  

##### STEP1 - QC  
**command:** 
```
lcWGSpipeline_step1-fastqc.py -p EXAMPLE.ckpt
```  
**flag(s):** -p or --ckpt_file  

**description:** The next step takes the checkpoint file created in the last step and prepares scripts to quality-check the raw fastqs. I recommended that you download the multiQC .html file to your local machine and open it in an internet brower to check for overall raw data quality. **Sedna's installation of multiQC is a mamba package. See Sedna's docs for setting up multiQC on your account.**

**output:** 
1) `./scripts/EXAMPLE-raw_fastqcARRAY.sh` (this file requires an array input file `./scripts/EXAMPLE-raw_fqcARRAY_input.txt`)
2) `./scripts/EXAMPLE-raw_multiqcSLURM.sh`.

**success message:** Step 1 has finished successfully! You will find two new scripts in `./scripts/`: `EXAMPLE-raw_fastqcARRAY.sh` and `EXAMPLE-raw_multiqcSLURM.sh`. `EXAMPLE-raw_fastqcARRAY.sh` must run prior to launching `EXAMPLE-raw_multiqcSLURM.sh`. However, there is no need to wait for these jobs to finish to move to step 2. To continue on, call step 2 to generate a script for TRIMMOMATIC. Remember to pass the checkpoint (.ckpt) file with the '-p' flag."  

**results:** 
1) a report for each fastq in `./fastqc/raw/`
2) multiqc_report.html  

##### STEP2 - TRIM  
**command:** 
```
lcWGSpipeline_step2-trim.py -p EXAMPLE.ckpt [--disable_fastp]
```  
**flag(s):** -p or --ckpt_file; --disable_fastp (set this flag if you would like to skip fastp (the tool for removing polyG tails)  

**description:** The next step takes the checkpoint file created in the last step and prepares scripts to trim adapters from raw fastqs, optionally clip polyG tails, and quality-check the trimmed/clipped fastqs. I recommended that you download the multiQC .html file to your local machine and open it in an internet brower to check for overall data quality before continuing to alignment.  **this will be updated in later versions of the pipeline, but for now once your scripts are created you need to go into the `EXAMPLE-trim_multiqcSLURM.sh` script and update the path to where your download of multiqc is. It currently is set up for just Laura's directory**

**output:**
1) `./scripts/EXAMPLE_trimARRAY.sh` (with array input file `./scripts/EXAMPLE_trimARRAY_input.txt`)
2) `./scripts/EXAMPLE-trim_fastqcARRAY.sh` (with array input file `./scripts/EXAMPLE-trim_fqcARRAY_input.txt`)
3) `./scripts/EXAMPLE-trim_multiqcSLURM.sh`.  

**success message:** "Step 2 has finished successfully! You will find three new scripts in `./scripts/`: `EXAMPLE_trimARRAY.sh`, `EXAMPLE-trim_fastqcARRAY.sh`, `EXAMPLE-trim_multiqcSLURM.sh`. Each script must run in the above order and each job must finish before submitting the next. While there is no need to wait before running step 3 to generate the alignment script, it is probably wise to wait until the multiQC script has completed and the results have been viewed in a web browser prior to submitting the script written by step 3. Remember to pass the checkpoint (.ckpt) file with the '-p' flag."  

**results:** 
`./trimmed/` will contain
1) a trimmed, paired file for each fastq
2) a trimmed, unparied file for each fastq
3) a timmed_clipped_paired file for each fastq (if --disable_fastp was **not** set)
4) a report for each fastq in `./fastqc/trimmed/`
5) multiqc_report_1.html  


##### STEP3 - ALIGN  
**command:** 
```
lcWGSpipeline_step3-align.py -p EXAMPLE.ckpt
```

**flag(s):** -p or --ckpt_file

**description:** This step aligns reads to a reference genome with BWA and runs the aligned reads through SAMTOOLS: 'fixmate' cleans up the read pairings and flags from BWA; a pair of 'view' statements converts the .sam file to a .bam file and filters the .bam file for non-unique and poor quality mappings; and 'sort' sorts the read pairings by coordinate (instead of read name). After a .bam file is built, duplicate reads are removed and (if the data is PE) overlapping reads are clipped to generate the final .bam. This step concludes by creating a file of depths for each alignment and generating a new file (`EXAMPLE_bamslist.txt`), which lists all final .bam files that will serve as input to ANGSD. In order for the `depthsARRAY.sh` script to run, make sure that `mean_cov_ind.py` is in the same bin directory as your lcWGSpipeline scripts.

**output:** 
1) `./scripts/EXAMPLE_alignARRAY.sh` (with array input file `./scripts/EXAMPLE_alignARRAY_input.txt`)
2) `./scripts/EXAMPLE_depthsARRAY.sh` (with array input file `./scripts/EXAMPLE_depthsARRAY_input.txt`)
3) `./EXAMPLE_bamslist.txt`

**success message:** Step 3 has finished successfully! You will find two new scripts in `./scripts/`: `EXAMPLE_alignARRAY.sh` and `EXAMPLE_depthsARRAY.sh`. These scripts must be run in the order given above and both must finish before moving on to step 4. After `EXAMPLE_depthsARRAY.sh` has run, download `bamtools/EXAMPLE_depths.csv` and generate a barchart of mean depths by individual. This will help determine whether any individuals fall substantially below the average depth (usually, we use a cutoff of 1x). If you identify samples with coverage that is 'too low', add the sample id(s) to a new file, referred to as the 'blacklist' of individuals to be excluded from genotype likelhiood calculations and the final data sets. After generating this exclude-list, you can continue to step 4 to write scripts for generating the final data sets. Remember to pass the checkpoint (.ckpt) file with the '-p' flag AND the exclude-list file with the '-b' flag.

**results:** 
1) a sorted, deduplicated, and clipped BAM file for each sample in ./bamtools
2) a .bai file for each sample in ./bamtools
3) ./bamtools/EXAMPLE_depths.tsv

##### STEP4 - RAW-GLS  
**command:** 
```
lcWGSpipeline_step4-rawGLS.py \
    -p EXAMPLE.ckpt \
    -e excludelist_1x.txt \
    -q 15 \
    -i 10 \
    -d 1 \
    -D 5 \
    -x q15d1D5 \
    --global_gls
```  

**flag(s):** 
--ckpt_file or -p: the checkpoint file created in step 0.
--exclude_inds or -e: a file listing any individuals that should be removed from downstream analyses
--global_gls: set this flag if you would like to calculate genotype likelihoods for all sites across the genome (if this flag is not set, gls will only be calculated for SNPs)
--quality_val or -q: value for minQ and minMapQ (default = 15)
--minInd or -i: use this flag to specify a minInd parameter (not required)
--minDepth_factor or -d: a factor to define minDepth: minDepth_factor * N = minDepth (default = 1)
--maxDepth_factor or -D: a factor to define maxDepth: maxDepth_factor * N = maxDepth (default = 5)
--prefix_extension or -x: option to add to the prefix to distinguish between parameterizations (not required)

**description:** This step calculates genotype likelihoods for putatively polymorphic sites (and all sites, if --global_gls is set), parallelized by chromosome, and concatenates results.

**output:** 
1) `./scripts/EXAMPLE_polymorphicARRAY.sh` (with array input file `./scripts/EXAMPLE_angsdARRAY_input.txt`)
2) a file of the included (not excluded) bams (`EXAMPLE_filtered_bamslist.txt`)
3) `./scripts/EXAMPLE_concatenate_polymorphic_beagles.sh`
4) `./scripts/EXAMPLE_concatenate_polymorphic_mafs.sh`
If you specified --global_gls:
5) `./scripts/EXAMPLE_globalARRAY.sh`
6) `./scripts/EXAMPLE_concatenate_global_beagles.sh`
7) `./scripts/EXAMPLE_concatenate_global_mafs.sh`
If you specified a prefix extension:
8) `./EXAMPLE-ADDITION.ckpt`

**success message:** Step 4 has finished successfully!
`EXAMPLE_polymorphicARRAY.sh` calculates genotype likelihoods across all polymorphic sites on each chromosome (separately).")
After this has run, you will have genotype likelihoods (gls) and allele frequencies (maf) for all putatively variable sites (SNPs).
`EXAMPLE_concatenate_polymorphic_beagles.sh` concatenates beagles from all SNPs across all chromosomes.
`EXAMPLE_concatenate_polymorphic_mafs.sh` concatenates mafs from all SNPs across all chromosomes, generating and indexing a sites file for the whole genome.
`EXAMPLE_polymorphicARRAY.sh` must have finished prior to running concatenation, but concatenation scripts can run in parallel.")
(If you specified --global_gls:) `EXAMPLE_globalARRAY.sh` calculates genotype likelihoods across all sites on each chromosome (separately).")
Both `EXAMPLE_polymorphicARRAY.sh` and `EXAMPLE_globalARRAY.sh` can run simultaneously.
After they have run, you will have genotype likelihoods (gls) and allele frequencies (maf) for all sites in the genome (global) and putatively variable sites (polymorphic).
`EXAMPLE_concatenate_global_beagles.sh` concatenates beagles from all positions across all chromosomes.
`EXAMPLE_concatenate_global_mafs.sh` concatenates mafs from all positions across all chromosomes, generating and indexing a sites file for the whole genome.
`EXAMPLE_globalARRAY.sh` must have finished prior to running concatenation, but concatenation scripts can run in parallel.
(If you specified a prefix extension:) You have specified a prefix extension. To facilitate downstream branching, a new checkpoint file has been written:
	`EXAMPLE-ADDITION.ckpt` contains all the information from EXAMPLE.ckpt
Additionally, the parameterization used for calculating genotype likelihoods in ANGSD has been recorded in `EXAMPLE-ADDITION.ckpt`. Including:  
	minDepth <default 1>  
	maxDepth <default 5>  
	minQ <default 15>  
(If the you set minInd:)  
	minInd <value>  

**results**  
1\) `./gls/EXAMPLE-ADDITION_{chrom}_polymorphic.beagle.gz` \(one for each chromosome\)  
2\) `./gls/EXAMPLE-ADDITION_{chrom}_polymorphic.mafs.gz` \(one for each chromosome\)  
3\) `./gls/EXAMPLE-ADDITION_wholegenome-polymorphic.beagle.gz`  
4\) `./gls/EXAMPLE-ADDITION_wholegenome-polymorphic.mafs.gz`  
5\) `./gls/EXAMPLE-ADDITION_wholegenome-polymorphic.sites`  
6\) `./gls/EXAMPLE-ADDITION_wholegenome-polymorphic.sites.bin`  
7\) `./gls/EXAMPLE-ADDITION_wholegenome-polymorphic.sites.idx`  
If you specified --global_gls:  
8\) `./gls/EXAMPLE-ADDITION_{chrom}_global.beagle.gz` \(one for each chromosome\)  
9\) `./gls/EXAMPLE-ADDITION_{chrom}_global.mafs.gz` \(one for each chromosome\)  
10\) `./gls/EXAMPLE-ADDITION_wholegenome-global.beagle.gz`  
11\) `./gls/EXAMPLE-ADDITION_wholegenome-global.mafs.gz`  
12\) `./gls/EXAMPLE-ADDITION_wholegenome-global.sites`  
13\) `./gls/EXAMPLE-ADDITION_wholegenome-global.sites.bin`  
14\) `./gls/EXAMPLE-ADDITION_wholegenome-global.sites.idx`  

##### STEP5 - FILTERED-GLS
**command:** 
```
lcWGSpipeline_step5-filteredGLS.py -p EXAMPLE-ADDITION.ckpt
```  
**flag(s):** 
--ckpt_file or -p: the checkpoint file created in step 0.
--alpha or -a: significance threshold for Chi-Square tests
--parallelize_n or -n: a number of parallel jobs to array (this value must be <999; if a value >999 is provided, the value will be overwritten with 999)
--ngsParalog_executable or -x: the full path to the ngsParalog executable
--sigTest or -s: the full path to ngsParalog_sigTest.R (available on github within this project)

**description:** This step begins by dividing the rawGLS sites file into a number of subsets (these files are named "subsitesa", "subsitesb"...). After creating these subsets, samtools mpileup runs and output is fed directly to ngsParalog to calculate likelihood ratios of a site being paralogous (mpileup files are not printed out). The resulting lr files are concatenated, resulting in one file of likelihood ratio values for the entire genome. The R script, ngsParalog_sigTest.R, runs a chisq test for each site to determine whether it is statistically likely to be paralogous (applying a Bonferroni correction). The R script returns a file of sites to retain for analysis, which is indexed and passed to ANGSD to, functionally, re-run step4 exclusively over the sites that survived filtering (putativey homologous sites).

**output:**
1) `./paralog/subsites{letter}` (letter = a, b, c,...parallelize_n)
2) `./scripts/EXAMPLE-ADDITION_paralogARRAY.sh`
2) `./scripts/EXAMPLE-ADDITION_sigLR.sh`
3) `./scripts/EXAMPLE-ADDITION_polymorphic-filteredARRAY.sh`
4) `./scripts/EXAMPLE-ADDITION_concatenate_wgph_beagles.sh`
5) `./scripts/EXAMPLE-ADDITION_concatenate_wgph_mafs.sh`

**success message:** Step 5 has finished successfully.
Note that 'wgph' abbreviates 'whole genome polymorphic homologous', representing files of paralog-filtered (putatively homologous) SNPs across the whole genome.
You will find a number of subset sites files, depending on the value you passed to parallelize_n.
There are five new scripts:
	`EXAMPLE-ADDITION_paralogARRAY.sh` generates an mpileup file for every subset sites file and runs ngsParalog.
	`EXAMPLE-ADDITION_sigLR.sh` concatenates the ngsParalog results from the subset sites files, identifies sites that are statistically likely to be homologous, and indexes those sites with angsd.
	`EXAMPLE-ADDITION_polymorphic-filteredARRAY.sh` is very similar to step4, calculating genotype likelihoods for homologous sites, parallelized by chromosome.
	`EXAMPLE-ADDITION_concatenate_wgph_beagles.sh` concatenates the beagle files resulting from executing `EXAMPLE-ADDITION_polymorphic-filteredARRAY.sh`
	`EXAMPLE-ADDITION_concatenate_wgph_mafs.sh` concatenates the mafs files resulting from executing `EXAMPLE-ADDITION_polymorphic-filteredARRAY.sh`
Both concatenation scripts can run simultaneously.
Note that a whole-genome sites file is not generated from the whole-genome mafs file (as in step4), because the whole-genome sites file results from `EXAMPLE-ADDITION_sigLR.sh`

**results:**  
1) ./paralog/subsites{letter}.lr (letter = a, b, c,...parallelize_n)
2) ./paralog/EXAMPLE-ADDITION_wholegenome.lr
3) ./paralog/EXAMPLE-ADDITION_wholegenome_retain.sites
4) ./paralog/EXAMPLE-ADDITION_wholegenome_retain.sites.bin
5) ./paralog/EXAMPLE-ADDITION_wholegenome_retain.sites.idx
6) ./gls/EXAMPLE-ADDITION_{chrom}_polymorphic-filtered.beagle.gz (one for each chromosome)
7) ./gls/EXAMPLE-ADDITION_{chrom}_polymorphic-filtered.mafs.gz (one for each chromosome)
8) ./gls/EXAMPLE-ADDITION_wgph.beagle.gz
9) ./gls/EXAMPLE-ADDITION_wgph.mafs.gz


### DATA ANALYSIS
##### PCA & ADMIXTURE  
**command:** `lcWGSpipeline_pca-admixture_vX.Y.py -p `EXAMPLE.ckpt` -k {10}`  

**flag(s):** -p or --ckpt_file; -k or --k_val_max (default = 10)  

**description:** Running this step results in a script to generate a new beagle file (a concatenation of all genotype likelihoods across all chromosomes), which will serve as input for scripts to execute PCA and ADMIXTURE analysis. This step also generates a script to run PCA over each chromosome.  

**output:** `EXAMPLE_concatenate_beagles.sh`, `EXAMPLE_pcangsdARRAY.sh` (which uses array input `EXAMPLE_angsdARRAY_input.txt`), and `EXAMPLE_wholegenome_pcangsd-admixture.sh`.  

**success message:** "Three scripts have been generated: `EXAMPLE_concatenate_beagles.sh`, `EXAMPLE_pcangsdARRAY.sh`, and `EXAMPLE_wholegenome_pcangsd-admixture.sh`. `EXAMPLE_concatenate_beagles.sh` must run first. After it has finished, both `EXAMPLE_pcangsdARRAY.sh` and `EXAMPLE_wholegenome_pcangsd-admixture.sh` can run simultaneously. `EXAMPLE_pcangsdARRAY.sh` will run pca for each chromosome. `EXAMPLE_wholegenome_pcangsd-admixture.sh` will execute pca and admixture analysis for the whole genome polymorphic data."  

**results:**  The pcangsd scripts will return a covariance matrix that does not include labels for which sample corresponds to which row. The output is in the same order as the bam file order that went into the likelihood calculation. Check the "filtered_bamslist.txt" file for the order. 


### POPULATION-LEVEL ANALYSES  
**command:** `lcWGSpipeline_pop-analyses_vX.Y.py -p {EXAMPLE.ckpt} -g {pops.txt}`  
**flag(s):** -p or --ckpt_file; -g or --group_file (formatted as individual_ID<TAB>population)  
**description:** This step generates a script to calculate a number of population-pair metrics. Site allele frequencies are calculated for each population, serving as input for the other calculations. Metrics include FST, dxy, and expected heterozygosity.  
**output:** {EXAMPLE_popARRAY.sh} (which uses array input {EXAMPLE_angsdARRAY_input.txt}).  
**success message:** "A new script has been generated to calculate population-level and populaiton-pair-level metrics: {EXAMPLE_popARRAY.sh}."  
**results:**  

### HMM  
**command:** `lcWGSpipeline_hmm_vX.Y.py -p {EXAMPLE.ckpt} -r {/home/ltimm/bin/HMM_log10FST+1_3norm.R}`  
**flag(s):** -p or --ckpt_file; -r or --r_script_loc (full path to the R script HMM_log10FST+1_3norm.R)  
**description:** This step generates an array script which uses a Hidden Markov Model (HMM) to classify sites by fst value (high, low, background). This step runs over the whole genome for each population pair.  
**output:** {EXAMPLE_hmmARRAY.sh} (which uses array input {EXAMPLE_hmmARRAY_input.txt}).  
**success message:** "A new script has been generated to run HMM in R. Submitting {EXAMPLE_hmmARRAY.sh} will classify each site as high (1), low (2), or background differentiation (3), based on the fst value.  
Note: for this step to run correctly, you must have HMM_log10FST+1_3norm.R available for the script to access (in your path)."  
**results:**  

### LOSTRUCT  
**command:** `lcWGSpipeline_lostruct_vX.Y.py -p {EXAMPLE.ckpt} -w {100000} -n {3} -r {/home/ltimm/bin/<put-the-r-script-location-here}`  
**flag(s):** -p or --ckpt_file; -w or --window_size (the desired window size to run pcangsd over); -n or --num_pcs (number of principal components to retain for each window); -r or --rscript_loc (full path to the R script that converts cov to eigenvalues; if you grab the R script from the same repo you found this, this refers to "cov_to_eigen.R")  
**description:** Because lostruct does not accept genotype likelihood data as input, this step takes each chromosome of GL data (the polymorphic beagle files associated with the project) and divides them into windows of user-specified length. PCAngsd runs over each window and passes the output to the cov_to_eigen.R script (from Nicolas Lou) along with the window size and the number of principal components to retain.  
**output:** {EXAMPLE_eigenARRAY.sh} (which uses array input {EXAMPLE_eigenARRAY_input.txt}).  
**success message:** This is forthcoming...  
**results:**  
**Note:** for this step to run correctly, you must have cov_to_eigen.R available for the script to access (in your path).  

### ADAPTIVE REGIONS  
**command:**  
**flag(s):**  
**description:**  
**output:**  
**success message:**  
**results:**  

### LINKAGE  
**command:**  
**flag(s):**  
**description:**  
**output:**  
**success message:**  
**results:**  

### GENOTYPE HEATMAPS  
**command:**  
**flag(s):**  
**description:**  
**output:**  
**success message:**  
**results:**  
