#0) Introductions

#snakemake pipeline for QC of vcf files. Created first on 191115.

#1) Set up Variables from Config files
shell.prefix('source activate 191115_vcf_QC_snakemake_env;')

#Read in File Reading options from config file
inFolder = config['inFolder']
outFolder = config['outFolder']
tempFolder = config['tempFolder']
statsFolder = config['statsFolder']
gzvcf_suffix = config['gzvcf_suffix']

#Read in Options from config file
splitFinalVCF = config['splitFinalVCF']
blacklist = config['blacklist']
blacklist_file = config['blacklist_file']

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

#Declare alleles variable to establish a later used wildcard
alleles = ['Biallelic']
#alleles = ['Biallelic','Triallelic'] #old alleles variable used when also considering

#allFiles finds all gzvcf files in the specific input directory. Parses names. For example: input/name.vcf.gz becomes name
allFiles = [x.split(inFolder) for x in glob.glob(inFolder + "/*" + gzvcf_suffix)]
allFiles = [x[1] for x in allFiles]
allFiles = [x.split(gzvcf_suffix) for x in allFiles]
allFiles = [x[0] for x in allFiles]

#Finds files in allFiles that have a specified chromosome (excluding chrX and chrY). These are valid files for the QC pipeline
regex = "_chr(?:[1,2][1,2,3,4,5,6,7,8,9,0]|[1,2,3,4,5,6,7,8,9])\." #regex is a regular expression to find strings such as "_chr1_", "_chr17_", etc.
chr = [re.findall(regex,x) for x in allFiles] #finds specific chromosome for each file. If there is no line matching regex, then returns nothing
valid = np.where([len(x) > 0 for x in chr])[0].tolist() #returns positions in allFiles where there are valid files for QC
allFiles_valid = [allFiles[x] for x in valid] #allFiles valid only takes the files in directory that A) have the gzvcf suffic, and B) have a "_chr#_" string within them. This includes chromosomes 1-22, but excludes X and Y

chr = [re.findall(regex,x) for x in allFiles_valid] #retakes the chr file, but excluding files with chrX and chrY, as well as excluding other files in the input directory not fulfilling criteria

#samples1 is the sample name before "_chr#_" and samples2 is the samples name after
allFiles_split = [re.split(regex,x) for x in allFiles_valid]
samples1 = [x[0] for x in allFiles_split]
samples2 = [x[1] for x in allFiles_split]

#convert samples1 and samples2 to sets to only keep unique values
samples1_set = set(samples1)
samples2_set = set(samples2)

#find indices for all unique values in both samples1 and samples2
samples1_indices = [[i for i, x in enumerate(samples1) if x == y] for y in samples1_set]
samples2_indices = [[i for i, x in enumerate(samples2) if x == y] for y in samples2_set]

#takes all intersections between samples1_indices and samples2_indices to return coordinates for every unique pair
def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3
combos = [intersection(x,y) for x,y in list(iter.product(samples1_indices,samples2_indices))]
combos = [x for x in combos if x != []]

#Makes an easy to access list for wildcards with unique sample IDs and the chromosomes associated with each
#at the end of this samples is in the following format: [[prefix1,postfix1,chrList1],[prefix2,postfix2,chrList2],...]
cohorts = []
for i in range(len(combos)):
    temp = []
    temp.append(samples1[combos[i][0]])
    temp.append(samples2[combos[i][0]])
    temp.append([chr[j] for j in combos[i]])
    cohorts.append(temp)

#Performs QC for each cohort in the input folder
for cohort in range(len(cohorts)):

    #Edits chromosome list to allow wildcards to work more easily
    chrs = [x[0] for x in cohorts[cohort][2]]
    chrs = [x[1:-1] for x in chrs]

    #3) "rule all" is a pseudorule that tells snakemake to make this output for each input
    #different final output for if you specify that the final vcf should be split or not. Splitting the vcf will invoke an extra rule.
    if splitFinalVCF:
        rule all:
            input:
                expand(outFolder + '{sample1}_{chr}.{sample2}_{allele}_QCFinished.recode.vcf.gz', allele = alleles, sample1 = cohorts[cohort][0], sample2 = cohorts[cohort][1], chr = chrs)
    else:
        rule all:
            input:
                expand(outFolder + '{sample1}_chrAll.{sample2}_{allele}_QCFinished.recode.vcf.gz', allele = alleles, sample1 = cohorts[cohort][0], sample2 = cohorts[cohort][1])

    #4) Filter based on Blacklist
    rule blacklist:
        input:
            vcfgz = inFolder + '{sample1}_{chr}.{sample2}' + gzvcf_suffix,
            blacklist_file = blacklist_file
        output:
            blacklist_filtered = tempFolder + '{sample1}_{chr}.{sample2}_blacklistFiltered' + gzvcf_suffix,
            stats = statsFolder + '{sample1}_{chr}.{sample2}_stats.txt'
        shell:
            "module load vcftools/0.1.15;"
            "module load bcftools/1.9;"
            "bcftools stats {input.vcfgz} > {output.stats};"
            "vcftools --gzvcf {input.vcfgz} --exclude-positions {input.blacklist} --stdout --recode --recode-INFO-all | bgzip -c > {output.blacklist_filtered};"

    #5) Chunk. Duplicated the rule based on whether theres a blacklist step or not because of different inputs. Can compact, but placing a function within input to change the list is unwieldly
    if Filter_Blacklist:
        rule chunk:
            input:
                vcfgz = tempFolder + '{sample1}_{chr}.{sample2}_blacklistFiltered' + gzvcf_suffix,
                chrLengths = chromosomeLengths
            output:
                chunked = tempFolder + '{sample1}_{chr}.{sample2}_{chunk}' + gzvcf_suffix,
                stats_input = statsFolder + '{sample1}_{chr}.{sample2}_blacklistFiltered_stats.txt',
                stats_output = statsFolder + '{sample1}_{chr}.{sample2}_{chunk}_stats.txt',
            params:
                chunkString = '', #empty string. Commands below fills this so that it may be fed to the shell command
                tempFolder = tempFolder,
                regex = regex,
            run:
                #Extracts the chromosome value from the output folder
                filename = output.chunked.split(params.tempFolder)
                filename = filename[1]
                chromosome = re.findall(params.regex,filename)[0]
                chromosome = chromosome[1:-1]

                #Extracts the length of the chromosome for the output file chromosome
                with open(input.chrLengths) as sizes:
                    contents=list(csv.reader(sizes,delimiter='\t'))
                contents = list(map(list, zip(*contents)))
                chrLen = int(contents[1][contents[0].index(chromosome)])

                #Extracts the chunk # from the output filename
                chunkNum = output.chunked.split("_chunk")[1]
                chunkNum = chunkNum.split(gzvcf_suffix)[0]

                #calculates the genetic loci for this chunk from the chromosome size and chunk number
                loci1 = (int(chunkNum)-1)*math.ceil(chrLen/NUM_CHUNK)+1
                loci2 = int(chunkNum)*math.ceil(chrLen/NUM_CHUNK)

                #Uses tabix to split vcf file to extract only this chunks portion of the chromosome
                params.chunkString = chromosome+':'+ str(loci1) +'-'+ str(loci2)
                shell("module load bcftools/1.9;tabix -p vcf {input.vcfgz};tabix -h {input.vcfgz} {params.chunkString} > {output.chunked};")

                #performs stats
                shell("module load bcftools/1.9;bcftools stats {input.vcfgz} > {output.stats_input};bcftools stats {output.chunked} > {output.stats_output};")
    else:
        rule chunk:
            input:
                vcfgz = inFolder + '{sample1}_{chr}.{sample2}' + gzvcf_suffix,
                chrLengths = chromosomeLengths
            output:
                chunked = tempFolder + '{sample1}_{chr}.{sample2}_{chunk}' + gzvcf_suffix,
                stats_input = statsFolder + '{sample1}_{chr}.{sample2}_stats.txt',
                stats_output = statsFolder + '{sample1}_{chr}.{sample2}_{chunk}_stats.txt',
            params:
                chunkString = '', #empty string. Commands below fills this so that it may be fed to the shell command
                tempFolder = tempFolder,
                regex = regex,
            run:
                #Extracts the chromosome value from the output folder
                filename = output.chunked.split(params.tempFolder)
                filename = filename[1]
                chromosome = re.findall(params.regex,filename)[0]
                chromosome = chromosome[1:-1]

                #Extracts the length of the chromosome for the output file chromosome
                with open(input.chrLengths) as sizes:
                    contents=list(csv.reader(sizes,delimiter='\t'))
                contents = list(map(list, zip(*contents)))
                chrLen = int(contents[1][contents[0].index(chromosome)])

                #Extracts the chunk # from the output filename
                chunkNum = output.chunked.split("_chunk")[1]
                chunkNum = chunkNum.split(gzvcf_suffix)[0]

                #calculates the genetic loci for this chunk from the chromosome size and chunk number
                loci1 = (int(chunkNum)-1)*math.ceil(chrLen/NUM_CHUNK)+1
                loci2 = int(chunkNum)*math.ceil(chrLen/NUM_CHUNK)

                #Uses tabix to split vcf file to extract only this chunks portion of the chromosome
                params.chunkString = chromosome+':'+ str(loci1) +'-'+ str(loci2)
                shell("module load bcftools/1.9;tabix -p vcf {input.vcfgz};tabix -h {input.vcfgz} {params.chunkString} > {output.chunked};")

                #performs stats
                shell("module load bcftools/1.9;bcftools stats {input.vcfgz} > {output.stats_input};bcftools stats {output.chunked} > {output.stats_output};")

    #6) Separate biallelics and triallelics for later on
    rule Filter0_separateBiallelics:
        input:
            chunked = tempFolder + '{sample1}_{chr}.{sample2}_{chunk}' + gzvcf_suffix
        output:
            biallelic_Filter0 = tempFolder + '{sample1}_{chr}.{sample2}' + '_Biallelic_{chunk}_Filter0.recode.vcf.gz',
            stats = statsFolder + '{sample1}_{chr}.{sample2}' + '_Biallelic_{chunk}_Filter0_stats.txt'
        shell:
            "module load vcftools/0.1.15;"
            "module load bcftools/1.9;"
            "vcftools --gzvcf {input.chunked} --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --stdout | bgzip -c > {output.biallelic_Filter0};"
            #"bcftools view --threads 5 -m2 -M2 -v snps -Oz {input.vcfgz} > {output.biallelic_Filter0};" #may speed up QC pipeline slightly. Haven't tested
            "bcftools stats {output.biallelic_Filter0} > {output.stats};"

    #Removed triallelic steps. Stretch Goal
    # rule Filter0_separateTriallelics:
    #     input:
    #         vcfgz = inFolder + '{sample1}_{chr}.{sample2}' + gzvcf_suffix
    #     output:
    #         triallelic_Filter0 = outFolder + '{sample1}_chrAll_{sample2}_Triallelic_QCFinished.recode.vcf',
    #         stats = statsFolder + '{sample1}_{chr}.{sample2}' + '_Triallelic_Filter0_stats.txt'
    #     params:
    #         unsplit_file = tempFolder + '{sample1}_{chr}.{sample2}' + '_Triallelic_unsplit.recode.vcf.gz'
    #     shell:
    #         "module load vcftools/0.1.15; module load bcftools/1.9;"
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
            Filter0 = tempFolder + '{sample1}_{chr}.{sample2}_{allele}_{chunk}_Filter0.recode.vcf.gz'
        output:
            Filter1 = tempFolder + '{sample1}_{chr}.{sample2}_{allele}_{chunk}_Filter1.recode.vcf.gz',
            stats = statsFolder + '{sample1}_{chr}.{sample2}_{allele}_{chunk}_Filter1_stats.txt'
        shell:
            "module load vcftools/0.1.15;"
            "module load bcftools/1.9;"
            "vcftools --gzvcf {input.Filter0} --remove-filtered-all --stdout --recode --recode-INFO-all | bgzip -c > {output.Filter1};"
            #"bcftools view --threads 5 -f PASS {input.Filter0} > {output.Filter1};" #may speed up QC pipeline slightly. Haven't tested
            "bcftools stats {output.Filter1} > {output.stats};"

    #8) Filter for genotype level depth (DP) (GDP)
    rule Filter2_GDP:
        input:
            Filter1 = tempFolder + '{sample1}_{chr}.{sample2}_{allele}_{chunk}_Filter1.recode.vcf.gz'
        output:
            Filter2 = tempFolder + '{sample1}_{chr}.{sample2}_{allele}_{chunk}_Filter2.recode.vcf.gz',
            stats = statsFolder + '{sample1}_{chr}.{sample2}_{allele}_{chunk}_Filter2_stats.txt'
        params:
            GDP_thresh = GDP
        shell:
            "module load vcflib/v1.0.0-rc0; module load bcftools/1.9;"

            #tabix must be created in this step and not step before because of how snakemake temporally creates files. Because tabix does not have an "output", it would be created before the output from the previous rule if placed within the rule.
            "tabix -p vcf {input.Filter1};"
            "vcffilter -g \"DP > {params.GDP_thresh}\" {input.Filter1} | bgzip -c > {output.Filter2};"
            "bcftools stats {output.Filter2} > {output.stats};"


    #9) Filter for Genome Quality (GQ)
    rule Filter3_GQ:
        input:
            Filter2 = tempFolder + '{sample1}_{chr}.{sample2}_{allele}_{chunk}_Filter2.recode.vcf.gz'
        output:
            Filter3 = tempFolder + '{sample1}_{chr}.{sample2}_{allele}_{chunk}_Filter3.recode.vcf.gz',
            stats = statsFolder + '{sample1}_{chr}.{sample2}_{allele}_{chunk}_Filter3_stats.txt'
        params:
            GQ_thresh = GQ
        shell:
            "module load vcflib/v1.0.0-rc0; module load bcftools/1.9;"
            "tabix -p vcf {input.Filter2};"
            "vcffilter -g \"GQ > {params.GQ_thresh}\" {input.Filter2} | bgzip -c > {output.Filter3};"
            "bcftools stats {output.Filter3} > {output.stats};"

    #10) Filter for SNP missingness
    rule Filter4_SNP_Missingess:
        input:
            Filter3 = tempFolder + '{sample1}_{chr}.{sample2}_{allele}_{chunk}_Filter3.recode.vcf.gz'
        output:
            Filter4 = tempFolder + '{sample1}_{chr}.{sample2}_{allele}_{chunk}_Filter4.recode.vcf.gz',
            stats = statsFolder + '{sample1}_{chr}.{sample2}_{allele}_{chunk}_Filter4_stats.txt'
        params:
            MISS_THRESH_SNP = MISS_THRESH_SNP
        shell:
            "module load vcftools/0.1.15; module load bcftools/1.9;"
            "vcftools --gzvcf {input.Filter3} --max-missing {params.MISS_THRESH_SNP} --stdout --recode --recode-INFO-all | bgzip -c > {output.Filter4};"
            "bcftools stats {output.Filter4} > {output.stats};"

    #11) Filter for Overall Read Depth (DP) (ODP)
    rule Filter5_ODP:
        input:
            Filter4 = tempFolder + '{sample1}_{chr}.{sample2}_{allele}_{chunk}_Filter4.recode.vcf.gz'
        output:
            Filter5 = tempFolder + '{sample1}_{chr}.{sample2}_{allele}_{chunk}_Filter5.recode.vcf.gz',
            stats = statsFolder + '{sample1}_{chr}.{sample2}_{allele}_{chunk}_Filter5_stats.txt'
        params:
            ODP_thresh = ODP
        shell:
            "module load vcflib/v1.0.0-rc0; module load bcftools/1.9;"
            "tabix -p vcf {input.Filter4};"
            "vcffilter -f \"DP > {params.ODP_thresh}\" {input.Filter4} | bgzip -c > {output.Filter5};"
            "bcftools stats {output.Filter5} > {output.stats};"

    #12) Filter for Mapping Quality (MQ)
    rule Filter6_MQ:
        input:
            Filter5 = tempFolder + '{sample1}_{chr}.{sample2}_{allele}_{chunk}_Filter5.recode.vcf.gz'
        output:
            Filter6 = tempFolder + '{sample1}_{chr}.{sample2}_{allele}_{chunk}_Filter6.recode.vcf.gz',
            stats = statsFolder + '{sample1}_{chr}.{sample2}_{allele}_{chunk}_Filter6_stats.txt'
        params:
            MQ_MAX_THRESH = MQ_MAX,
            MQ_MIN_THRESH = MQ_MIN
        shell:
            "module load vcflib/v1.0.0-rc0; module load bcftools/1.9;"
            "tabix -p vcf {input.Filter5};"
            "vcffilter -f \"MQ > {params.MQ_MIN_THRESH} & MQ < {params.MQ_MAX_THRESH}\" {input.Filter5} | bgzip -c > {output.Filter6};"
            "bcftools stats {output.Filter6} > {output.stats};"

    #13) Separate out Indel files and SNP files in Biallelic files.
    rule Biallelic_Separate_Indels_and_SNPs:
        input:
            Filter6 = tempFolder + '{sample1}_{chr}.{sample2}_Biallelic_{chunk}_Filter6.recode.vcf.gz'
        output:
            Filter6_SNPs = tempFolder + '{sample1}_{chr}.{sample2}_Biallelic_{chunk}_Filter6_SNPs.recode.vcf.gz',
            Filter6_SNPs_stats = statsFolder + '{sample1}_{chr}.{sample2}_Biallelic_{chunk}_Filter6_SNPs_stats.txt',
            Filter6_Indels = tempFolder + '{sample1}_{chr}.{sample2}_Biallelic_{chunk}_Filter6_Indels.recode.vcf.gz',
            Filter6_Indels_stats = statsFolder + '{sample1}_{chr}.{sample2}_Biallelic_{chunk}_Filter6_Indels_stats.txt'
        shell:
            "module load vcftools/0.1.15; module load bcftools/1.9;"

            "vcftools --gzvcf {input.Filter6} --remove-indels --stdout --recode --recode-INFO-all | bgzip -c > {output.Filter6_SNPs};"
            "vcftools --gzvcf {input.Filter6} --keep-only-indels --stdout --recode --recode-INFO-all | bgzip -c > {output.Filter6_Indels};"

            "bcftools stats {output.Filter6_SNPs} > {output.Filter6_SNPs_stats};"
            "bcftools stats {output.Filter6_Indels} > {output.Filter6_Indels_stats};"


    #14) Filter Biallelic SNPs for VQSLOD
    rule Biallelic_SNPs_Filter7_VQSLOD:
        input:
            Filter6_SNPs = tempFolder + '{sample1}_{chr}.{sample2}_Biallelic_{chunk}_Filter6_SNPs.recode.vcf.gz'
        output:
            Filter7_SNPs = tempFolder + '{sample1}_{chr}.{sample2}_Biallelic_{chunk}_Filter7_SNPs.recode.vcf.gz',
            stats = statsFolder + '{sample1}_{chr}.{sample2}_Biallelic_{chunk}_Filter7_SNPs_stats.txt'
        params:
            VQSLOD_thresh = VQSLOD
        shell:
            "module load vcflib/v1.0.0-rc0; module load bcftools/1.9;"
            "tabix -p vcf {input.Filter6_SNPs};"
            "vcffilter -f \"VQSLOD > {params.VQSLOD_thresh}\" {input.Filter6_SNPs} | bgzip -c > {output.Filter7_SNPs};"
            "bcftools stats {output.Filter7_SNPs} > {output.stats};"

    #15) Combine Biallelic Indels and SNPs
    rule Biallelic_Combine_Indels_and_SNPs:
        input:
            Filter7_SNPs = tempFolder + '{sample1}_{chr}.{sample2}_Biallelic_{chunk}_Filter7_SNPs.recode.vcf.gz',
            Filter6_Indels = tempFolder + '{sample1}_{chr}.{sample2}_Biallelic_{chunk}_Filter6_Indels.recode.vcf.gz'
        output:
            Filter7 = tempFolder + '{sample1}_{chr}.{sample2}_Biallelic_{chunk}_Filter7.recode.vcf.gz',
            stats = statsFolder + '{sample1}_{chr}.{sample2}_Biallelic_{chunk}_Filter7_stats.txt'
        shell:
            "module load bcftools/1.9;"
            "bcftools concat {input.Filter6_Indels} {input.Filter7_SNPs} | bcftools sort | bgzip -c > {output.Filter7};"
            "bcftools stats {output.Filter7} > {output.stats};"

    #16) Filter via Inbreeding_Coef
    rule Filter8_Inbreeding_Coef:
        input:
            Filter7 = tempFolder + '{sample1}_{chr}.{sample2}_{allele}_{chunk}_Filter7.recode.vcf.gz'
        output:
            Filter8 = tempFolder + '{sample1}_{chr}.{sample2}_{allele}_{chunk}_Filter8.recode.vcf.gz',
            stats = statsFolder + '{sample1}_{chr}.{sample2}_{allele}_{chunk}_Filter8_stats.txt'
        params:
            INBREEDING_COEF = INBREEDING_COEF
        run:
            shell("module load vcflib/v1.0.0-rc0; module load bcftools/1.9;tabix -p vcf {input.Filter7};")

            #Negative values for inbreeding coeficient are expressed as (0 - [value]) due to how vcffilter is built. Hence, an if statement to separate negative from positive input INBREEDING_COEF values
            if params.INBREEDING_COEF >= 0:
                shell("module load vcflib/v1.0.0-rc0; module load bcftools/1.9;vcffilter -f \"InbreedingCoeff > {params.INBREEDING_COEF}\" {input.Filter7} |bgzip -c > {output.Filter8};")
            else:
                params.INBREEDING_COEF = -params.INBREEDING_COEF
                shell("module load vcflib/v1.0.0-rc0; module load bcftools/1.9;vcffilter -f \"InbreedingCoeff > ( 0 - {params.INBREEDING_COEF} )\" {input.Filter7} |bgzip -c > {output.Filter8};")

            shell("module load bcftools/1.9;bcftools stats {output.Filter8} > {output.stats};")

    #17) recombines chunks and chromosomes in one step
    rule recombine_Chromosomes:
        input:
            separate = expand(tempFolder + '{sample1}_{chr}.{sample2}_{allele}_{chunk}_Filter8.recode.vcf.gz', chunk = chunks,
            allele = alleles, sample1 = cohorts[cohort][0], sample2 = cohorts[cohort][1], chr = chrs)
        output:
            recombined = tempFolder + '{sample1}_chrAll_{sample2}_{allele}_Filter8.recode.vcf.gz',
            stats = statsFolder + '{sample1}_chrAll_{sample2}_{allele}_Filter8_stats.txt'
        params:
            concat_string = '',
        run:
            for i in input.separate:
                params.concat_string = params.concat_string+str(i)+' '

            shell("module load vcftools/0.1.15; module load bcftools/1.9;vcf-concat {params.concat_string} | bgzip > {output.recombined};")
            shell("module load bcftools/1.9; bcftools stats {output.recombined} > {output.stats};")

    rule sort:
        input:
            Filter8 = tempFolder + '{sample1}_chrAll_{sample2}_{allele}_Filter8.recode.vcf.gz'
        output:
            Filter8_sort = tempFolder + '{sample1}_chrAll_{sample2}_{allele}_Filter8_sorted.recode.vcf.gz',
            stats = statsFolder + '{sample1}_chrAll_{sample2}_{allele}_Filter8_sorted_stats.txt'
        params:
            memory = cluster.mem
        shell:
            "module load bcftools/1.9;
            "bcftools sort {input.Filter8} -Oz -m {params.memory} > {output.Filter8_sort};"
            "bcftools stats {output.Filter8_sort} > {output.stats};"

    #18) Filter via Sample level Missingness
    rule Filter9_Sample_Missingness:
        input:
            Filter8 = tempFolder + '{sample1}_chrAll_{sample2}_{allele}_Filter8_sorted.recode.vcf.gz'
        output:
            Filter9 = tempFolder + '{sample1}_chrAll_{sample2}_{allele}_Filter9.recode.vcf.gz',
            stats = statsFolder + '{sample1}_chrAll_{sample2}_{allele}_Filter9_stats.txt'
        params:
            MISS_THRESH_INDI = MISS_THRESH_INDI,
            missingness_name = tempFolder + '{sample1}_chrAll_{sample2}_{allele}_Filter8',
            missingness_filtered = tempFolder + '{sample1}_chrAll_{sample2}_{allele}_Filter8_imiss_filtered.txt'
        shell:
            "module load vcftools/0.1.15; module load bcftools/1.9;module load R/3.5.3;"

            #create sample missingness file
            "vcftools --gzvcf {input.Filter8} --missing-indv --out {params.missingness_name};"

            #filter the sample missingness file for significance
            "Rscript Accessory_Scripts/filterMissingnessIndividualFile.R {params.missingness_name}.imiss {params.MISS_THRESH_INDI} {params.missingness_filtered};"

            #filter the samples with significant missingness out of vcf
            "vcftools --gzvcf {input.Filter8} --remove {params.missingness_filtered} --stdout --recode --recode-INFO-all | bgzip -c > {output.Filter9};"

            "bcftools stats {output.Filter9} > {output.stats};"

    #19) Filter via Relatedness
    rule Filter10_Relatedness:
        input:
            Filter9 = tempFolder + '{sample1}_chrAll_{sample2}_{allele}_Filter9.recode.vcf.gz'
        output:
            Final = outFolder + '{sample1}_chrAll.{sample2}_{allele}_QCFinished.recode.vcf.gz',
            stats = statsFolder + '{sample1}_chrAll_{sample2}_{allele}_QCFinished_stats.txt'
        params:
            RELATEDNESS_THRESH = RELATEDNESS_THRESH,
            relatedness_name = tempFolder + '{sample1}_chrAll_{sample2}_{allele}_Filter9',
            relatedness_filtered = tempFolder + '{sample1}_chrAll_{sample2}_{allele}_Filter9_relatedness2_filtered.txt'
        shell:
            "module load vcftools/0.1.15; module load bcftools/1.9;module load R/3.5.3;"

            #create relatedness file
            "vcftools --gzvcf {input.Filter9} --relatedness2 --out {params.relatedness_name};"

            #filter the relatedness file for significance
            "Rscript Accessory_Scripts/filterRelatednessFile.R {params.relatedness_name}.relatedness2 {params.RELATEDNESS_THRESH} {params.relatedness_filtered};"

            ##filter the samples with significant relatedness out of vcf. The RScript optimizes number of samples retained. More details about algorithm inside Rscript
            "vcftools --gzvcf {input.Filter9} --remove {params.relatedness_filtered} --stdout --recode --recode-INFO-all | bgzip -c > {output.Final};"

            "bcftools stats {output.Final} > {output.stats};"

    #20) Splits the chromosomes again if specified in config.
    rule Split_ChrAll:
        input:
            Final = outFolder + '{sample1}_chrAll.{sample2}_{allele}_QCFinished.recode.vcf.gz'
        output:
            Final_split = outFolder + '{sample1}_{chr}.{sample2}_{allele}_QCFinished.recode.vcf.gz'
        params:
            chromosome = '',
            regex = regex,
            outFolder = outFolder
        run:
            filename = output.Final_split.split(params.outFolder)
            filename = filename[1]

            chromosome2 = re.findall(params.regex,filename)[0]
            params.chromosome = chromosome2[1:-1]
            shell("module load bcftools/1.9;tabix -p vcf {input.Final};tabix {input.Final} {params.chromosome}> {output.Final_split};")
