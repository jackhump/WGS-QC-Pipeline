# WGS QC Pipeline

A snakemake pipeline for performing QC on WGS VCF files. Includes QCs for Minor Allele Frequency (MAF), Hardy Weinberg Equilibrium (HWE), Relatedness, and Missingness (sample and SNP). Takes in a gzipped vcf file.

Based on script in Adelson et al. 2019

Currently version: 1.2

# Command

source snakejob [path to config.yaml]

# Config Options


General Options

*splitFinalVCF*: If true, splits the final output vcf.gz to create per chromosome vcf.gz files

*blacklist*: If true, filters initially based on blacklist file provided

*blactlist_file*: Path to blacklist file. hg38 blacklist is contained within Accessory_Scripts folder


Chunking Metrics

*NUM_CHUNK*: How much to break up the chromosome. For example, if 4, every chromosome will be broken in four and processed in parallel. Too high values may cause bugs

*chromosomeLengths*: File of chromosome lengths. Used for chunking




File Configs

*inFolder*: Folder where input files are located

*outFolder*: Folder where final output is stored

*interFolder*: Folder where intermediate files are stored. These are files such as relatedness, missingness, and other metrics.

*statsFolder*: Folder where stats are stored.

*gzvcf_suffix*: Suffix for gzipped vcf file. Used to extract file names.




QC Metrics

*MISS_THRESH_SNP*: Threshold of missingness of a SNP, above which we exlude SNP from analysis. Default is 0.15

*ODP*: Overall Read Depth. Default is 5000. Paper's default is 25000. Should be dependent on total number of samples in cohort, and coverage. 

*MQ_MIN*: Mapping Quality lower threshold. Default is 58.75

*MQ_MAX*: Mapping Quality upper threshold. Default is 61.25

*VQSLOD*: Filtering criteria applied only to SNVs. Default is 7.81

*GDP*: Genotype level read depth. Default is 10.

*GQ*: Genome Quality. Default is 20.

*Inbreeding Coefficient*: Default is -0.8

*MISS_THRESH_INDI*: Threshold of missingness of a Sample, above which we exclude Sample from analysis. Default is .1

*RELATEDNESS_THRESH*: Threshold of relatedness, above which we exclude a Sample based on relation to other Samples. Default is .125

# To do

Enable option to make intermediate outputs temporary
Include plotting options
Shift from vcf-concat to bcftools concat (might be faster)
Chunk bug fix (referenced below)
If blacklist is true, but no path to a file is provided, use default

# Known bugs

High chunk number will cause empty files to be written. At some point, shift chunking over to SnpSift to chunk based on file size rather than chromosome location. Will also speed up chunking related processes. 

# Conda recipe

conda create –n 191115_vcf_QC_snakemake_env

conda config --add channels bioconda

conda config --add channels defaults

conda config --add channels conda-forge

conda install snakemake

# Filter Levels

Filter-1: Filter for provided blacklist, if any

Chunk: Chunks chromosomes, if chunking value is set to >1

Filter0: Separate out Biallelics

Filter1: Filter for GATK Passes

Filter2: Filter for genotype level depth (DP) (GDP)

Filter3: Filter for Genome Quality (GQ)

Filter4: Filter for SNP missingness

Filter5: Filter for Overall Read Depth (DP) (ODP)

Filter6: Filter for Mapping Quality (MQ)

Filter7: Filter SNPs for VQSLOD

Filter8: Filter via Inbreeding_Coef

Filter9: Filter via Sample level Missingness

Filter10: Filter via Relatedness