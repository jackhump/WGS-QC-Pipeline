# WGS QC Pipeline

A snakemake pipeline for performing QC on WGS VCF files. Includes QCs for Minor Allele Frequency (MAF), Hardy Weinberg Equilibrium (HWE), Relatedness, and Missingness (sample and SNP). Takes in a gzipped vcf file.

# Command

source snakejob [path to config.yaml]

# Config Options

File Configs
inFolder: Folder where input files are located
outFolder: Folder where final output is stored
interFolder: Folder where intermediate files are stores. These are files such as relatedness, missingness, and other metrics.
gzvcf_suffix: Suffix for gzipped vcf file. Used to extract file names.

QC Metrics
MISS_THRESH_SNP: Threshold of missingness of a SNP, above which we exlude SNP from analysis. Default is 0.15
MISS_THRESH_INDI: Threshold of missingness of a Sample, above which we exclude Sample from analysis. Default is .1
RELATEDNESS_THRESH: Threshold of relatedness, above which we exclude a Sample based on relation to other Samples. Default is .125
MAF: Minor Allele frequency threshold. We keep all Minor Allele Frequencies greater than this. Default is .01
HWE: Hardy Weinberg Equilibrium threshold. Default is .00000001

# Conda recipe

conda create –n 191115_vcf_QC_snakemake_env
conda config --add channels bioconda
conda config --add channels defaults
conda config --add channels conda-forge
conda install snakemake
