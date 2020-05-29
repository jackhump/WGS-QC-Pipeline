# WGS QC Pipeline
------------------
Written by Jack Humphrey and Sanan Venkatesh

A snakemake pipeline for performing quality control (QC) on whole genome sequencing (WGS) data. 
Includes filters for Genotype quality (MQ, DP, GQ, VQSLOD), minor allele frequency, relatedness, and missingness (sample and variant). 

Heavily inspired by *Adelson, R.P., Renton, A.E., Li, W. et al. Empirical design of a variant quality control pipeline for whole genome sequencing data using replicate discordance. Sci Rep 9, 16156 (2019)*

## Dependencies and Installation

This pipeline is set to work with the Mount Sinai HPC 'Minerva' which uses the LSF queuing system and the Lmod module environment system. It is assumed that the following dependencies are present:

- python 3.6.8 (modules snakemake, yaml, pandas, numpy, itertools, math, csv, re, os)
- R 3.6.0 (packages dplyr, readr, stringr)
- bcftools 1.9
- king/2.1.6
- vcftools/0.1.15
- vcflib/v1.0.0-rc0
- plink2

All of these dependencies are present on Minerva except for snakemake. Install this through conda:

```
ml anaconda3/latest
conda create –n WGS-QC-pipeline
conda config --add channels bioconda
conda config --add channels defaults
conda config --add channels conda-forge
conda install snakemake
conda activate WGS-QC-pipeline
```


## Inputs

A set of jointly called VCFs from GATK, split by chromosome. Each file must have "chr1","chr2", etc in the file name.

The full paths to each file should be in a text file. 

Variants *must* have been called using GATK. Check whether your data has VQSLOD, DP, MQ and InbreedingCoef fields.

All config options are stored in `config.yaml`.

dbSNP RS IDs are provided by ensembl. To generate the annotation VCF, run: 

```
cd data/; sh create_ensembl_vcf.sh
```


## Outputs

Using the `dataCode` parameter specified in the config.yaml, the following subfolders will be created:

- <dataCode>/input/
- <dataCode>/temp/
- <dataCode>/output/
- <dataCode>/relatedness/
- <dataCode>/stats/
 
The initial VCF files will be symlinked to `input/`. Any intermediary files will be written to `temp/`. Outputs from KING on the relatedness between samples will be written to `relatedness/`. Statistics on numbers of variants from each step will be written to `stats/`. Final output files will be written to `output`.

The output files will be:

- **chrAll_QCFinished_full.anno** - the full set of QC'd variants, in VCF and PLINK format, annotated with dbSNP RS IDs
- **chrAll_QCFinished_MAF<threshold>.anno** - the minor allele frequency filtered set of QC'd variants, in VCF and PLINK format, annotated with dbSNP RS IDs
- **all_variant_stats_collated.txt** - the total number of variants retained at each step of the pipeline.
 

## Specifying samples to remove

VCF IDs for samples to be removed from the analysis should be in a text file. The relative path to this file should be specified in the config.yaml under the key `removeSamples`



# Parallel execution on MSSM HPC

```
./snakejob -c [path to config.yaml] <-n> <-a> 
```

`-n` specifies dry run mode. Snakemake will list the steps that it plans to execute.
`-a` specifies the account code for running jobs on the cluster


# Config Options

kept inside config.yaml

### General Options

**dataCode**: _string_

The name of your analysis. All files will be written to subfolders of a folder with this name.

**vcfFileList**: _FILE_

Path to file containing full paths to each joint called VCF file.

**samples_to_remove**: _FILE_, False

Path to file containing VCF sample IDs of any samples to be removed. Otherwise set as False.

**blacklist_file**: _FILE_ 

Path to blacklist file. hg38 blacklist is contained within `/data` folder. 
Current blacklist taken from https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg38-blacklist.v2.bed.gz; ref: https://www.nature.com/articles/s41598-019-45839-z

**liftOverhg19hg38**: _True,False_

whether to lift over the variants from hg19 to hg38. This requires the human genome fasta file `hg38.fa` to be present in `data/`. Symlink this from wherever you keep a copy.

**NUM_CHUNK**: _number_ 

How many pieces to break up the chromosome. This allows for parallel execution. Chunk numbers of 4 have been tested and work without error. Larger numbers may cause errors due to the creation of empty chunks 

**chromosomeLengths**: _FILE_ 

Text file of chromosome lengths. Used for chunking. hg38 chromosome sizes is contained within `data` folder. Sizes were taken from IGV at https://github.com/igvteam/igv/blob/master/genomes/sizes/hg38.chrom.sizes

**splitFinalVCF**: _True,False_ 

If true, splits the final output vcf.gz to create per chromosome vcf.gz files - CURRENTLY HARDCODED OFF

### QC Thresholds

**MISS_THRESH_SNP**: _number_ 

Threshold of missingness of a SNP, above which we exlude SNP from analysis. Default is 0.15

**ODP**: _number_

Overall Read Depth. Default is 5000. Number is dependent on total number of samples in cohort multiplied by the coverage. See below for more details.

**MQ_MIN**: _number_

Mapping Quality lower threshold. Default is 58.75

**MQ_MAX**: _number_

Mapping Quality upper threshold. Default is 61.25

**VQSLOD**: _number_ 

Variant quality log-odds score from GATK. Applied only to SNPs. Default is 7.81

**GDP**: _number_

Genotype level read depth. Default is 10.

**GQ**: _number_

Genotype Quality. Default is 20.

**Inbreeding Coefficient**: _number_

Used instead of Hardy-Weinberg. Default is -0.8.

**MISS_THRESH_INDI**: _number_

Individual missingness threshold. Default i 0.1.

**RELATEDNESS_THRESH**: _number_

Threshold of relatedness, above which we exclude a Sample based on relation to other Samples. Default is .125


# How to estimate DP, MQ and VQSLOD thresholds from your own data

using one of your autosomal VCF files (say chr22 for example), create diagnostic plots for the VCF parameters MQ (mapping quality), DP (total read depth) and VQSLOD (variant quality score log-odds):

```
cd scripts/
sh get_quality_metrics <vcf.gz> <data code>
```

This will create the folder scripts/<data code> which will contain the raw numbers for each metric along with histograms of each metric.

In our experience, MQ and VQSLOD do not change between datasets but DP does as it corresponds to the total read depth at each base - the average DP being the average per-sample read depth (usually 30X) multiplied by the sample number. 

So 300 samples at 30x coverage should have a median DP of around 9000.


# Known bugs

High chunk number will cause empty files to be written. At some point, shift chunking over to SnpSift to chunk based on file size rather than chromosome location. Will also speed up chunking related processes.

Please report any bugs as issues on this github.


# Filter Levels

Filter-1: Filter for provided blacklist and any known samples for removal.

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

Filter11: Filter by Minor Allele Frequency
