#0) Introductions

#snakemake pipeline for QC of vcf files. Created first on 191115.
# Sanan Venkatesh and Jack Humphrey
print( " * WGS QC PIPELINE * ")
print( " authors: Jack Humphrey and Sanan Venkatesh")
print( " inspired by Adelson et al, 2019")

#1) Set up Variables from Config files
R_VERSION = "R/4.0.3"
# should be set to some default conda environment that contains snakemake
#shell.prefix('ml anaconda3/latest; conda activate WGS-QC-pipeline;')

#Read in File Reading options from config file
dataCode = config['dataCode']

inFolder = dataCode + "/input/"
outFolder = dataCode + "/output/"
tempFolder = dataCode + "/temp/"
statsFolder = dataCode + "/stats/"
relatedFolder = dataCode + "/relatedness/"

#Read in Options from config file
splitFinalVCF = config['splitFinalVCF']
#blacklistFile = config['blacklistFile']
removeSamples = config['removeSamples']
liftOverhg19hg38 = config['liftOverhg19hg38']

genomeBuild = config['genome_build']
print(" * Genome build used is: ", genomeBuild)
if genomeBuild == "hg38":
    genomeFASTA = "data/hg38.fa"
    blacklistFile = "data/hg38-blacklist.v2.bed.gz"

if genomeBuild == "hg19":
    genomeFASTA = "data/hg19.fa"
    blacklistFile = "data/hg19-blacklist.v2.bed.gz"

print(" * using blacklist:" , blacklistFile)
#Read in Chunking Metrics from config file
NUM_CHUNK = config['NUM_CHUNK']
chunks = ['chunk'+str(x) for x in range(1,NUM_CHUNK+1)]
chromosomeLengths = config['chromosomeLengths']

#Read in QC Metrics from config file
MISS_THRESH_SNP = config['MISS_THRESH_SNP']
ODP = config['ODP']
MQ_MAX = config['MQ_MAX']
MQ_MIN = config['MQ_MIN']
VQSLOD = config['VQSLOD']
GDP = config['GDP']
GQ = config['GQ']
INBREEDING_COEF = config['INBREEDING_COEF']
MISS_THRESH_INDI = config['MISS_THRESH_INDI']
RELATEDNESS_THRESH = config['RELATEDNESS_THRESH']

#Makes sure there is no distinction whether you enter /input or /input/
if (inFolder.endswith("/") == False):
    inFolder = inFolder + "/"
if (outFolder.endswith("/") == False):
    outFolder = outFolder + "/"
if (tempFolder.endswith("/") == False):
    tempFolder = tempFolder + "/"
if (statsFolder.endswith("/") == False):
    statsFolder = statsFolder + "/"

#2) Pull Samples names from inFolder and uses string manipulation to parse their names. This does handle multiple cohorts within the same input folder

#import libraries
import glob
import os
import re
import numpy as np
import itertools as iter
import math
import csv
import subprocess
#Declare alleles variable to establish a later used wildcard
alleles = ['Biallelic']
#alleles = ['Biallelic','Triallelic'] #old alleles variable used when also considering

#allFiles finds all gzvcf files in the specific input directory. Parses names. For example: input/name.vcf.gz becomes name

#Finds files in allFiles that have a specified chromosome (excluding chrX and chrY). These are valid files for the QC pipeline

# replace Sanan's complicated code with simply:
#   find all vcf.gz files in {dataCode}/input/
#   keep only autosomes (chr1-22)
#   refer throughout as {dataCode}{chr}.whatever
vcf_file_list = config['vcfFileList']

# specify whether chr1 or 1
chr_format = config["chr_format"]

# assume all files have a "chrN" in the file name
chr_string = "chr[0-9]+"

regex = re.compile(chr_string)

vcfFiles = [line.rstrip('\n') for line in open(vcf_file_list)]

# remove entries that don't end with ".vcf.gz"
vcfFiles = [ v for v in vcfFiles if v.endswith(".vcf.gz")]

print(" * Found %s gzipped VCF files" % len(vcfFiles) )

# remove entries that don't have chr names
vcfFiles = [ v for v in vcfFiles if regex.search(v) is not None ]

print(" * Found %s VCF files with chromosomes between 1 and 22" % len(vcfFiles) )

print(" * Chunking: %s chunks selected" % NUM_CHUNK)


# get out chrs
chrs = [regex.search(vcf).group() for vcf in vcfFiles if regex.search(vcf) is not None] 
# sort files by chr number 
chr_n = [ int(i.split("chr")[1]) for i in chrs]
chr_sorted = sorted(zip(chr_n,chrs, vcfFiles) )

# get the sorted values back out
vcfFiles = [i[2] for i in chr_sorted]
chrs = [i[1] for i in chr_sorted]

# create list of new output file names
symlinkedFiles = [ inFolder + chr + "_input.vcf.gz" for chr in chrs ] 

if splitFinalVCF:
    final_output = expand(outFolder + '{chr}_QCFinished.vcf.gz', allele = alleles,  chr = chrs)
else:
    print(" * Output options: don't split by chromosome at the end")
    final_output = expand(outFolder + 'chrAll_QCFinished.vcf.gz' )

MAF_threshold = str(config["MAF_threshold"])

rule all:
        input:
            #expand( inFolder + "{chr}_input.vcf.gz", chr = chrs)
            #expand(tempFolder + '{chr}_{chunk}_chunked.vcf.gz', chr = chrs, chunk = chunks) # rule: chunk
            #expand(tempFolder + '{chr}_{chunk}_filtered.vcf.gz', chr = chrs, chunk = chunks) # rule: filterRegionsAndSamples
            #expand(tempFolder + '{chr}_{chunk}_Biallelic.recode.vcf.gz', chr = chrs, chunk = chunks) # rule: Filter0_separateBiallelics
            #expand(tempFolder + '{chr}_{chunk}_Filter1.recode.vcf.gz', chr = chrs, chunk = chunks) # rule: Filter1_GATK_PASS
            #expand(tempFolder + '{chr}_{chunk}_Filter2.recode.vcf.gz', chr = chrs, chunk = chunks) # rule: Filter2_GDP
            #expand(tempFolder + '{chr}_{chunk}_Filter3.recode.vcf.gz', chr = chrs, chunk = chunks) # rule: Filter3_GQ
            #expand(tempFolder + '{chr}_{chunk}_Filter4.recode.vcf.gz', chr = chrs, chunk = chunks) # rule: Filter4_SNP_Missingess
            #expand(tempFolder + '{chr}_{chunk}_Filter5.recode.vcf.gz', chr = chrs, chunk = chunks) # rule: Filter5_ODP
            #expand(tempFolder + '{chr}_{chunk}_Filter6.recode.vcf.gz', chr = chrs, chunk = chunks) # rule: Filter6_MQ
            #expand(tempFolder + '{chr}_{chunk}_Filter6_SNPs.recode.vcf.gz', chr = chrs, chunk = chunks) # rule: Biallelic_Separate_Indels_and_SNPs
            #expand(tempFolder + '{chr}_{chunk}_Filter7_SNPs.recode.vcf.gz', chr = chrs, chunk = chunks) # rule: Biallelic_SNPs_Filter7_VQSLOD
            #expand(tempFolder + '{chr}_{chunk}_CombineSNPsIndels.recode.vcf.gz', chr = chrs, chunk = chunks) # rule: Biallelic_Combine_Indels_and_SNPs
            #expand(tempFolder + '{chr}_{chunk}_Filter8.recode.vcf.gz', chr = chrs, chunk = chunks) # rule: Filter8_Inbreeding_Coef
            #expand(tempFolder + '{chr}_Filter8.recode.vcf.gz', chr = chrs) # rule: recombineChunks
            #expand(tempFolder + 'chrAll_Recombined.vcf.gz') # rule: recombineChromosomes
            #expand(tempFolder + 'chrAll_Recombined.IDs.vcf.gz') # rule: setID
            #expand(tempFolder + 'chrAll_Recombined.bed') # rule: convertVCFtoPLINK
            #expand(tempFolder + 'chrAll_Filter9.bed') # rule: Filter9_Sample_Missingness
            #expand(relatedFolder + 'chrAll_QCFinishedunrelated.txt') # rule: KingRelatedness
            #expand(outFolder + 'chrAll_QCFinished_full.bed') # rule: removeRelatedSamples
            #expand(outFolder + 'chrAll_QCFinished_MAF' + MAF_threshold + ".bed") # rule: filterMAF
            #expand(outFolder + 'chrAll_QCFinished_{file}.vcf.gz', file = ['full', 'MAF' + MAF_threshold ]) # rule: convertPlinkToVCF
            #expand(outFolder + 'chrAll_QCFinished_{file}.anno.vcf.gz', file = ['full', 'MAF' + MAF_threshold ]) # rule: annotateVCF
            outFolder + "all_variant_stats_collated.txt" # rule: collateStats

# for each file in the VCF file - symlink to input folder
rule symlinkVCFs:
    output:
        expand( inFolder + "{chr}_input.vcf.gz", chr = chrs),
        expand( inFolder + "{chr}_input.vcf.gz.tbi", chr = chrs)
    run:
    # symlink files to named file to inFolder/newFiles[i]
        for i in range(len(vcfFiles)):
            os.symlink(os.path.abspath(vcfFiles[i]), symlinkedFiles[i], target_is_directory = False, dir_fd = None)
            os.symlink(os.path.abspath(vcfFiles[i] + ".tbi"), symlinkedFiles[i] + ".tbi", target_is_directory = False, dir_fd = None)

# CHUNKING

# chunk number set by user
# read in chrLengths
# for each chr and chunk - calculate
rule chunk:
    input:
        vcfgz = ancient(inFolder + "{chr}_input.vcf.gz"),
        tbi = ancient(inFolder + "{chr}_input.vcf.gz.tbi"),
        #vcfgz = tempFolder + '{chr}_filtered.vcf.gz',
        chrLengths = ancient(chromosomeLengths)
    output:
        chunked = temp(tempFolder + '{chr}_{chunk}_chunked.vcf.gz'),
        stats_output = statsFolder + '{chr}_{chunk}_Chunk_stats.txt'
    params:
        chunkString = ''
    message: "Creating chunk {wildcards.chunk} of {NUM_CHUNK} for {wildcards.chr}"
    run:
        import pandas as pd
        chromosome = wildcards.chr
        chunk = wildcards.chunk
        chunkNum = int(chunk.split("chunk")[1])

        chrSizes = pd.read_csv(chromosomeLengths, sep = "\t", index_col = 0, header = None)
        chrLen = chrSizes.loc[chromosome]

        if chr_format == "1":
            chromosome = chromosome.replace("chr", "")

        #calculates the genetic loci for this chunk from the chromosome size and chunk number
        start = (int(chunkNum)-1)*math.ceil(chrLen/NUM_CHUNK)+1
        end = int(chunkNum)*math.ceil(chrLen/NUM_CHUNK)

        #Uses tabix -f  to split vcf file to extract only this chunks portion of the chromosome
        params.chunkString = chromosome + ':' + str(start) + '-' + str(end)

        shell("ml bcftools/1.9; tabix -f -h {input.vcfgz} {params.chunkString} | bgzip -c > {output.chunked}; bcftools stats {output.chunked} > {output.stats_output};")

## STEP 2 - PER-CHUNK FILTERS

## LIFT OVER to HG38 if necessary
if liftOverhg19hg38 is True:
    print(" * VCF files will be lifted over from hg19 to hg38")
    filterRegionsSamplesInput = tempFolder + "{chr}_{chunk}_hg38_sorted.vcf.gz"
else:
    filterRegionsSamplesInput = tempFolder + "{chr}_{chunk}_chunked.vcf.gz"

# if VCFs are called with hg19
# lift over to hg38
rule liftOverVCFs:
    input:
        vcf = tempFolder + '{chr}_{chunk}_chunked.vcf.gz', 
        genome =  "data/hg38.fa"
    output:
        vcfgz_sorted = temp(tempFolder + "{chr}_{chunk}_hg38_sorted.vcf.gz")
    params:
        vcfgz = tempFolder + "{chr}_{chunk}_hg38.vcf.gz",
        memory = "45G",
        vcf =  tempFolder + "{chr}_{chunk}_hg38.vcf",
        chain = "data/hg19ToHg38.over.chain.gz"
    shell:
        "ml crossmap/0.3.2; ml bcftools/1.9;"
        "CrossMap.py vcf {params.chain} {input.vcf} {input.genome} {params.vcf} ;"
        "bgzip {params.vcf} ;"
        "bcftools sort {params.vcfgz} -Oz --max-mem {params.memory}  > {output.vcfgz_sorted};"
        "tabix {params.vcfgz};"
        #"tabix {output.vcfgz_sorted};"
        #"rm {params.vcfgz}"


if removeSamples is not False:
    print(" * Removing samples from file %s " % removeSamples)
    sample_filter_string = "--remove " + removeSamples
else:
    sample_filter_string = ""

# Filter SNPs based on Blacklist, Samples based on removeSamples list
rule filterRegionsAndSamples:
    input:
        vcfgz = ancient(filterRegionsSamplesInput),
        blacklist_file = blacklistFile
    output:
        blacklist_filtered = temp(tempFolder + '{chr}_{chunk}_filtered.vcf.gz'),
        stats1 = statsFolder + '{chr}_{chunk}_Initial_stats.txt',
        stats2 = statsFolder + '{chr}_{chunk}__BlacklistFiltered_stats.txt'
    params:
        sample_filter_string = sample_filter_string
    shell:
        "ml vcftools/0.1.15;"
        "ml bcftools/1.9;"
        "tabix -f {input.vcfgz};"
        "bcftools stats {input.vcfgz} > {output.stats1};"
        "vcftools --gzvcf {input.vcfgz} {params.sample_filter_string} --exclude-positions {input.blacklist_file} --stdout --recode --recode-INFO-all | bgzip -c > {output.blacklist_filtered};"
        "tabix -f  {output.blacklist_filtered};"
        "bcftools stats {output.blacklist_filtered} > {output.stats2};"

#6) Separate biallelics and triallelics for later on
rule Filter0_separateBiallelics:
    input:
        chunked = tempFolder + '{chr}_{chunk}_filtered.vcf.gz'
        #chunked = tempFolder + '{chr}_{chunk}_chunked.vcf.gz'
    output:
        biallelic_vcf = temp(tempFolder + '{chr}' + '_{chunk}_Biallelic.recode.vcf.gz'),
        stats = statsFolder + '{chr}' + '_{chunk}_separateBiallelic.stats.txt'
    #group: "per_chunk_filter"
    shell:
        "ml vcftools/0.1.15;"
        "ml bcftools/1.9;"
        "n_var=$(bcftools view -H {input.chunked} | wc -l);if [[ $n_var == 0 ]];then echo 'chunk is empty'; cp {input.chunked} {output.biallelic_vcf}; else "
        "vcftools --gzvcf {input.chunked} --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --stdout | bgzip -c > {output.biallelic_vcf};"
        "fi;"
        #"bcftools view --threads 5 -m2 -M2 -v snps -Oz {input.vcfgz} > {output.biallelic_vcf};" #may speed up QC pipeline slightly. Haven't tested
        "bcftools stats {output.biallelic_vcf} > {output.stats};"

#Removed triallelic steps. Stretch Goal
# rule Filter0_separateTriallelics:
#     input:
#         vcfgz = inFolder + '{chr}.vcf.gz'
#     output:
#         triallelic_vcf = outFolder + 'chrAll_Triallelic_QCFinished.recode.vcf',
#         stats = statsFolder + '{chr}' + '_Triallelic_Filter0_stats.txt'
#     params:
#         unsplit_file = tempFolder + '{chr}' + '_Triallelic_unsplit.recode.vcf.gz'
#     shell:
#         "ml vcftools/0.1.15; ml bcftools/1.9;"
#
#         #this would be where I say "--max-alleles 3" if we only want to keep triallelics
#         "vcftools --gzvcf {input.vcfgz} --min-alleles 3 --recode --recode-INFO-all --stdout | bgzip -c > {params.unsplit_file};"
#
#         #I might need to make it "-m -". I am unclear on what the doc means. People online just seem to do -m though, but theres a way to combine biallelics into multiallelics
#         #-m splits multiallelics into multiple lines of biallelics, -Oz means the output will be compressed vcf
#         "bcftools norm -m - -Oz {params.unsplit_file} -o {output.triallelic_Filter0};"
#         "bcftools stats {output.triallelic_Filter0} > {output.stats};"

#7) Keep sites that "PASS" by GATK
rule Filter1_GATK_PASS:
    input:
        vcf = tempFolder + '{chr}_{chunk}_Biallelic.recode.vcf.gz'
    output:
        vcf = temp(tempFolder + '{chr}_{chunk}_Filter1.recode.vcf.gz'),
        stats = statsFolder + '{chr}_{chunk}_Filter1_stats.txt'
    #group: "per_chunk_filter"
    shell:
        "ml vcftools/0.1.15;"
        "ml bcftools/1.9;"
        "nLine=$(bcftools view -H {input.vcf} | wc -l); if [ $nLine == 0 ]; then cp {input.vcf} {output.vcf}; touch {output.stats}; exit 0; fi;"
        "vcftools --gzvcf {input.vcf} --remove-filtered-all --stdout --recode --recode-INFO-all | bgzip -c > {output.vcf};"
        #"bcftools view --threads 5 -f PASS {input.vcf} > {output.vcf};" #may speed up QC pipeline slightly. Haven't tested
        "bcftools stats {output.vcf} > {output.stats};"

#8) Filter for genotype level depth (DP) (GDP)
rule Filter2_GDP:
    input:
        vcf = tempFolder + '{chr}_{chunk}_Filter1.recode.vcf.gz'
    output:
        vcf = temp(tempFolder + '{chr}_{chunk}_Filter2.recode.vcf.gz'),
        stats = statsFolder + '{chr}_{chunk}_Filter2_stats.txt'
    params:
        GDP_thresh = GDP
    #group: "per_chunk_filter"
    shell:
        "ml vcflib/v1.0.0-rc0; ml bcftools/1.9;"

        #tabix -f  must be created in this step and not step before because of how snakemake temporally creates files. Because tabix does not have an "output", it would be created before the output from the previous rule if placed within the rule.
        "tabix -f  -p vcf {input.vcf};"
        "nLine=$(bcftools view -H {input.vcf} | wc -l); if [ $nLine == 0 ]; then cp {input.vcf} {output.vcf}; touch {output.stats}; exit 0; fi;"
        "vcffilter -g \"DP > {params.GDP_thresh}\" {input.vcf} | bgzip -c > {output.vcf};"
        "bcftools stats {output.vcf} > {output.stats};"

#9) Filter for Genome Quality (GQ)
rule Filter3_GQ:
    input:
        vcf = tempFolder + '{chr}_{chunk}_Filter2.recode.vcf.gz'
    output:
        vcf = temp(tempFolder + '{chr}_{chunk}_Filter3.recode.vcf.gz'),
        stats = statsFolder + '{chr}_{chunk}_Filter3_stats.txt'
    params:
        GQ_thresh = GQ
    #group: "per_chunk_filter"
    shell:
        "ml vcflib/v1.0.0-rc0; ml bcftools/1.9;"
        "tabix -f  -p vcf {input.vcf};"
        "nLine=$(bcftools view -H {input.vcf} | wc -l); if [ $nLine == 0 ]; then cp {input.vcf} {output.vcf}; touch {output.stats}; exit 0; fi;"
        "vcffilter -g \"GQ > {params.GQ_thresh}\" {input.vcf} | bgzip -c > {output.vcf};"
        "bcftools stats {output.vcf} > {output.stats};"

#10) Filter for SNP missingness
rule Filter4_SNP_Missingess:
    input:
        vcf = tempFolder + '{chr}_{chunk}_Filter3.recode.vcf.gz'
    output:
        vcf = temp(tempFolder + '{chr}_{chunk}_Filter4.recode.vcf.gz'),
        stats = statsFolder + '{chr}_{chunk}_Filter4_stats.txt'
    params:
        MISS_THRESH_SNP = MISS_THRESH_SNP
    #group: "per_chunk_filter" 
    shell:
        "ml vcftools/0.1.15; ml bcftools/1.9;"
        "nLine=$(bcftools view -H {input.vcf} | wc -l); if [ $nLine == 0 ]; then cp {input.vcf} {output.vcf}; touch {output.stats}; exit 0; fi;"        
        "vcftools --gzvcf {input.vcf} --max-missing {params.MISS_THRESH_SNP} --stdout --recode --recode-INFO-all | bgzip -c > {output.vcf};"
        "bcftools stats {output.vcf} > {output.stats};"

#11) Filter for Overall Read Depth (DP) (ODP)
rule Filter5_ODP:
    input:
        vcf = tempFolder + '{chr}_{chunk}_Filter4.recode.vcf.gz'
    output:
        vcf = temp(tempFolder + '{chr}_{chunk}_Filter5.recode.vcf.gz'),
        stats = statsFolder + '{chr}_{chunk}_Filter5_stats.txt'
    params:
        ODP_thresh = ODP
    #group: "per_chunk_filter"
    shell:
        "ml vcflib/v1.0.0-rc0; ml bcftools/1.9;"
        "tabix -f  -p vcf {input.vcf};"
        "nLine=$(bcftools view -H {input.vcf} | wc -l); if [ $nLine == 0 ]; then cp {input.vcf} {output.vcf}; touch {output.stats}; exit 0; fi;"
        "vcffilter -f \"DP > {params.ODP_thresh}\" {input.vcf} | bgzip -c > {output.vcf};"
        "bcftools stats {output.vcf} > {output.stats};"

#12) Filter for Mapping Quality (MQ)
rule Filter6_MQ:
    input:
        vcf = tempFolder + '{chr}_{chunk}_Filter5.recode.vcf.gz'
    output:
        vcf = temp(tempFolder + '{chr}_{chunk}_Filter6.recode.vcf.gz'),
        stats = statsFolder + '{chr}_{chunk}_Filter6_stats.txt'
    params:
        MQ_MAX_THRESH = MQ_MAX,
        MQ_MIN_THRESH = MQ_MIN
    #group: "per_chunk_filter"
    shell:
        "ml vcflib/v1.0.0-rc0; ml bcftools/1.9;"
        "tabix -f  -p vcf {input.vcf};"
        "nLine=$(bcftools view -H {input.vcf} | wc -l); if [ $nLine == 0 ]; then cp {input.vcf} {output.vcf}; touch {output.stats}; exit 0; fi;" 
        #"vcffilter -f \"MQ > {params.MQ_MIN_THRESH} & MQ < {params.MQ_MAX_THRESH}\" {input.vcf} | bgzip -c > {output.vcf};"
        "bcftools filter -i \"MQ > {params.MQ_MIN_THRESH} && MQ < {params.MQ_MAX_THRESH}\" {input.vcf} | bgzip -c > {output.vcf};"
        "bcftools stats {output.vcf} > {output.stats};"

#13) Separate out Indel files and SNP files in Biallelic files.
rule Biallelic_Separate_Indels_and_SNPs:
    input:
        vcf = tempFolder + '{chr}_{chunk}_Filter6.recode.vcf.gz'
    output:
        vcf_SNPs = temp(tempFolder + '{chr}_{chunk}_Filter6_SNPs.recode.vcf.gz'),
        vcf_SNPs_stats = statsFolder + '{chr}_{chunk}_Filter6_SNPs_stats.txt',
        vcf_Indels = temp(tempFolder + '{chr}_{chunk}_Filter6_Indels.recode.vcf.gz'),
        vcf_Indels_stats = statsFolder + '{chr}_{chunk}_Filter6_Indels_stats.txt'
    #group: "per_chunk_filter"
    shell:
        "ml vcftools/0.1.15; ml bcftools/1.9;"
        "nLine=$(bcftools view -H {input.vcf} | wc -l); if [ $nLine == 0 ]; then cp {input.vcf} {output.vcf_SNPs}; cp {input.vcf} {output.vcf_Indels}; touch {output.vcf_SNPs_stats}; touch {output.vcf_Indels_stats}; exit 0; fi;"
        "vcftools --gzvcf {input.vcf} --remove-indels --stdout --recode --recode-INFO-all | bgzip -c > {output.vcf_SNPs};"
        "vcftools --gzvcf {input.vcf} --keep-only-indels --stdout --recode --recode-INFO-all | bgzip -c > {output.vcf_Indels};"

        "bcftools stats {output.vcf_SNPs} > {output.vcf_SNPs_stats};"
        "bcftools stats {output.vcf_Indels} > {output.vcf_Indels_stats};"

#14) Filter Biallelic SNPs for VQSLOD
rule Biallelic_SNPs_Filter7_VQSLOD:
    input:
        vcf = tempFolder + '{chr}_{chunk}_Filter6_SNPs.recode.vcf.gz'
    output:
        vcf = temp(tempFolder + '{chr}_{chunk}_Filter7_SNPs.recode.vcf.gz'),
        stats = statsFolder + '{chr}_{chunk}_Filter7_SNPs_stats.txt'
    params:
        VQSLOD_thresh = VQSLOD
    #group: "per_chunk_filter"
    shell:
        "ml vcflib/v1.0.0-rc0; ml bcftools/1.9;"
        "tabix -f -p vcf {input.vcf};"
        "nLine=$(bcftools view -H {input.vcf} | wc -l); if [ $nLine == 0 ]; then cp {input.vcf} {output.vcf}; touch {output.stats}; exit 0; fi;"
        "vcffilter -f \"VQSLOD > {params.VQSLOD_thresh}\" {input.vcf} | bgzip -c > {output.vcf};"
        "bcftools stats {output.vcf} > {output.stats};"

#15) Combine Biallelic Indels and SNPs
rule Biallelic_Combine_Indels_and_SNPs:
    input:
        SNPs = tempFolder + '{chr}_{chunk}_Filter7_SNPs.recode.vcf.gz',
        Indels = tempFolder + '{chr}_{chunk}_Filter6_Indels.recode.vcf.gz'
    output:
        vcf = temp(tempFolder + '{chr}_{chunk}_CombineSNPsIndels.recode.vcf.gz'),
        stats = statsFolder + '{chr}_{chunk}_CombineSNPsIndels_stats.txt'
    #group: "per_chunk_filter"
    shell:
        "ml bcftools/1.9;"
        "nLine=$(bcftools view -H {input.SNPs} | wc -l); if [ $nLine == 0 ]; then cp {input.SNPs} {output.vcf}; touch {output.stats}; exit 0; fi;"
        "bcftools concat {input.Indels} {input.SNPs} | bcftools sort | bgzip -c > {output.vcf};"
        "bcftools stats {output.vcf} > {output.stats};"

#16) Filter via Inbreeding_Coef
rule Filter8_Inbreeding_Coef:
    input:
        vcf = tempFolder + '{chr}_{chunk}_CombineSNPsIndels.recode.vcf.gz'
    output:
        vcf = temp(tempFolder + '{chr}_{chunk}_Filter8.recode.vcf.gz'),
        stats = statsFolder + '{chr}_{chunk}_Filter8_stats.txt'
    params:
        INBREEDING_COEF = INBREEDING_COEF
    #group: "per_chunk_filter"
    run:
        shell("ml vcflib/v1.0.0-rc0; ml bcftools/1.9;tabix -f -p vcf {input.vcf};")
        cmd = "ml bcftools/1.9; bcftools view -H " + input.vcf + "| wc -l"
        nLine = subprocess.check_output(cmd, shell = True)
        nLine = int(nLine.decode("utf-8").strip() )
        if nLine == 0:
            shell("cp {input.vcf} {output.vcf}; touch {output.stats}")
        else:     
            #Negative values for inbreeding coeficient are expressed as (0 - [value]) due to how vcffilter is built. Hence, an if statement to separate negative from positive input INBREEDING_COEF values
            if params.INBREEDING_COEF >= 0:
                shell("ml vcflib/v1.0.0-rc0; ml bcftools/1.9;vcffilter -f \"InbreedingCoeff > {params.INBREEDING_COEF}\" {input.vcf} |bgzip -c > {output.vcf};")
            else:
                params.INBREEDING_COEF = -params.INBREEDING_COEF
                shell("ml vcflib/v1.0.0-rc0; ml bcftools/1.9;vcffilter -f \"InbreedingCoeff > ( 0 - {params.INBREEDING_COEF} )\" {input.vcf} |bgzip -c > {output.vcf};")

            shell("ml bcftools/1.9;bcftools stats {output.vcf} > {output.stats};")

## STEP 3: AGGREGATE ALL CHUNKS TOGETHER

#17) recombines chunks and chromosomes in one step
# recombine chunks first and sort
# then merge chromosomes but do in numeric order - must be possible
rule recombineChunks:
    input:
        #expand("{sample}_{id}.txt", id=["1", "2"], allow_missing=True)
        expand(tempFolder +  "{chr}_{chunk}_Filter8.recode.vcf.gz", chunk = chunks, allow_missing=True) # expand only chunk and allele, not chr
    wildcard_constraints:
        chunk="chunk[0-9]+",
        chr="chr[0-9]+"
    output:
        temp(tempFolder + '{chr}_Filter8.recode.vcf.gz')
    shell:
        "ml vcftools/0.1.15; ml vcflib/v1.0.0-rc0; ml bcftools/1.9;"
        "bcftools concat {input} | bgzip -c > {output}"

rule recombineChromosomes:
    input:
        expand(tempFolder + '{chr}_Filter8.recode.vcf.gz', chr = chrs) # expand both chr and allele
    output:
        recombined = temp(tempFolder + 'chrAll_Recombined.vcf.gz'),
        stats = statsFolder + 'chrAll_Recombined_stats.txt'
    threads: 16
    shell:
        "ml bcftools/1.9;"
        "bcftools concat --threads {threads} -Oz -o {output.recombined} {input};"
        "bcftools stats {output.recombined} > {output.stats};"

# Set ID  to chr_pos - multiallelic SNPs have been removed already
# also make sure reference allele is right way around - currently deprecated
# using hg38.fa
# unless genome build is hg19
rule setID:
    input:
        vcf = tempFolder + 'chrAll_Recombined.vcf.gz',
        genome =  genomeFASTA
    output:
        vcf = tempFolder + 'chrAll_Recombined.IDs.vcf.gz'
    threads: 32
    shell:
        "ml bcftools/1.9;"
        #"bcftools norm -Ou --check-ref ws -f {input.genome} | "
        "bcftools annotate --threads {threads} -Oz --output {output.vcf} --set-id +'%CHROM\:%POS' {input.vcf}"

## STEP 4: QC ON ALL DATA TOGETHER IN PLINK

# make binary PLINK file
rule convertVCFtoPLINK:
    input:
        vcf = tempFolder + 'chrAll_Recombined.IDs.vcf.gz'
    output:
        bed = tempFolder + 'chrAll_Recombined.bed',
        bim =  tempFolder + 'chrAll_Recombined.bim',
        fam =  tempFolder + 'chrAll_Recombined.fam'
    params:
        prefix = tempFolder + 'chrAll_Recombined',
        #relatedness_filtered = tempFolder + 'chrAll_Filter9_relatedness2_filtered.txt'
    shell:
        "ml plink2;"
        "plink2 --vcf {input} --make-bed --out {params.prefix}  --output-chr chrM "
 
#18) Filter via Sample level Missingness
rule Filter9_Sample_Missingness:
    input:
        bed = tempFolder + 'chrAll_Recombined.bed',
        bim = tempFolder + 'chrAll_Recombined.bim',
        fam = tempFolder + 'chrAll_Recombined.fam'
        #vcf = tempFolder + 'all_Filter8_sorted.recode.vcf.gz'
    output:
        bed = tempFolder + 'chrAll_Filter9.bed',
        bim = tempFolder + 'chrAll_Filter9.bim',
        fam = tempFolder + 'chrAll_Filter9.fam'
    params:
        MISS_THRESH_INDI = MISS_THRESH_INDI,
        prefix = tempFolder + 'chrAll_Filter9',
    shell:
        "ml plink2;"
        "plink2  --bed {input.bed} --bim {input.bim} --fam {input.fam} --mind {MISS_THRESH_INDI} --out {params.prefix} --make-bed " # the following breaks KING --output-chr chrM "

# Calculate Relatedness
rule KingRelatedness:
    input:
        bed =  tempFolder + 'chrAll_Filter9.bed'
    output:
        samples_to_keep = relatedFolder + 'chrAll_QCFinishedunrelated.txt'
    params:
        prefix = relatedFolder + 'chrAll_QCFinished'
    shell:
        "ml king/2.1.6;"
        "king -b {input.bed} --duplicate --prefix {params.prefix};"
        "king -b {input.bed} --related --degree 3 --prefix {params.prefix};"
        "king -b {input.bed} --unrelated --degree 3 --prefix {params.prefix};"

# king outputs a set of individuals not related closer than third degree (a blood relative which includes the individual’s first-cousins, great-grandparents or great grandchildren)
# remove those individuals 
rule removeRelatedSamples:
    input:
        bed = tempFolder + 'chrAll_Filter9.bed',
        bim = tempFolder + 'chrAll_Filter9.bim',
        fam = tempFolder + 'chrAll_Filter9.fam',
        samples_to_keep = relatedFolder + 'chrAll_QCFinishedunrelated.txt'
    output:
        bed = outFolder + 'chrAll_QCFinished_full.bed',
        bim =  outFolder + 'chrAll_QCFinished_full.bim',
        fam =  outFolder + 'chrAll_QCFinished_full.fam'
    params:
        prefix =  outFolder + 'chrAll_QCFinished_full'
    shell:
        "ml plink2;"
        "plink2  --bed {input.bed} --bim {input.bim} --fam {input.fam} --keep {input.samples_to_keep} --out {params.prefix} --make-bed " #--output-chr chrM  "

# Filter on Minor Allele Frequency (MAF)
rule filterMAF:
    input:
       bed = outFolder + 'chrAll_QCFinished_full.bed',
       bim = outFolder + 'chrAll_QCFinished_full.bim',
       fam = outFolder + 'chrAll_QCFinished_full.fam'
    output:
       bed = outFolder + 'chrAll_QCFinished_MAF' + MAF_threshold + ".bed",
       bim = outFolder + 'chrAll_QCFinished_MAF' + MAF_threshold + ".bim",
       fam = outFolder + 'chrAll_QCFinished_MAF' + MAF_threshold + ".fam"
    params:
        prefix =  outFolder + 'chrAll_QCFinished_MAF' + MAF_threshold
    shell:
        "ml plink2;"
        "plink2 --bed {input.bed} --bim {input.bim} --fam {input.fam} --maf {MAF_threshold} --out {params.prefix} --make-bed " #--output-chr chrM  "

# Convert Full and MAF-filtered variant sets back to VCF
rule convertPlinkToVCF:
    input:
        bed = outFolder + 'chrAll_QCFinished_{file}.bed',
        bim = outFolder + 'chrAll_QCFinished_{file}.bim',
        fam = outFolder + 'chrAll_QCFinished_{file}.fam'
    output:
        vcf = outFolder + 'chrAll_QCFinished_{file}.vcf.gz',
        stats = statsFolder + 'chrAll_QCFinished_{file}_stats.txt'
    params:
        prefix =  outFolder + 'chrAll_QCFinished_{file}'
    shell:
        "ml plink2; ml bcftools/1.9;"
        "plink2 --bed {input.bed} --bim {input.bim} --fam {input.fam} --recode vcf bgz --out {params.prefix} --output-chr chrM  ;"
        "tabix {output.vcf};"
        "bcftools stats {output.vcf} > {output.stats}"

# add RS ID numbers to both VCFs
rule annotateVCF:
    input:
        vcf = outFolder + 'chrAll_QCFinished_{file}.vcf.gz'
    output:
        vcf = outFolder + 'chrAll_QCFinished_{file}.anno.vcf.gz'
    params:
        ensembl = "data/ensembl_v99_hg38.vcf.gz" 
    shell:
        "ml bcftools/1.9;"
        "bcftools annotate -Oz -o {output.vcf} -a {params.ensembl} -c ID {input.vcf};"
        "tabix {output.vcf}"


# put all stats outputs together to get numbers of variants at each stage
rule collateStats:
    input:
        expand(  outFolder + 'chrAll_QCFinished_{file}.vcf.gz', file = ['full', 'MAF' + MAF_threshold ] )
    output:
        outFolder + "all_variant_stats_collated.txt"
    params:
        script = "scripts/collate_all_stats.R"
    shell:
        "ml {R_VERSION};"
        "Rscript {params.script} {statsFolder} {output}" 
