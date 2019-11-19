#0) Introductions

#snakemake pipeline for QC of vcf files. Created first on 191115.

#1) Set up Variables from Config files
shell.prefix('export PS1="";source activate bazam-pipeline;')

inFolder = config['inFolder']
outFolder = config['outFolder']
interFolder = config['interFolder']

MISS_THRESH_SNP = config['MISS_THRESH_SNP']
MISS_THRESH_INDI = config['MISS_THRESH_INDI']
RELATEDNESS_THRESH = config['RELATEDNESS_THRESH']
MAF = config['MAF']
HWE = config['HWE']

gzvcf_suffix = config['gzvcf_suffix']

#2) Pull Samples names from inFolder #For now, only takes first file in folder. Debateably will expand to allow multiple vcf.gz inputs
import glob
import os

samples = [os.path.basename(x).strip(gzvcf_suffix) for x in glob.glob(inFolder + "/*" + gzvcf_suffix)]
#sample = sample[0]

#3) "rule all" is a pseudorule that tells snakemake to make this output for each input
rule all:
    input:
        expand(outFolder + '/{sample}_QCFinished.recode.vcf', sample = samples)
    #input: [outFolder + '/{sample}' + '_QCFinished.recode.vcf' for sample in samples]

#4) Recode the original vcf to filter for MAF, HWE. Also remove indels, and shrink alleles down to 2

rule Filter1:
    input:
        vcfgz = inFolder + '/{sample}' + gzvcf_suffix
    output:
        filter1 = interFolder + '/{sample}' + '_filter1.recode.vcf'
    params:
        output_name = interFolder + '/{sample}' + '_filter1'
    shell:
        "module load vcftools/0.1.15;"
        "vcftools --gzvcf {input.vcfgz} --out {params.output_name} --maf $MAF --hwe $HWE --remove-indels --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all;"

#5) Correct for Relatedness

rule Filter2:
    input:
        filter1 = interFolder + '/{sample}' + '_filter1.recode.vcf'
    output:
        filter2 = interFolder + '/{sample}' + '_filter2.recode.vcf'
    params:
        relatedness_file_name = interFolder + '/{sample}',
        relatedness_file = interFolder + '/{sample}' + '.relatedness2',
        relatedness_file_filter = interFolder + '/{sample}' + '_relatedness2_filtered.txt',
        output_name = interFolder + '/{sample}' + '_filter2',
        relatedness_thresh = RELATEDNESS_THRESH
    shell:
        "module load vcftools/0.1.15;"
        "vcftools --vcf {input.filter1} --out {params.relatedness_file_name} --relatedness2;"

        "./snakeRefScript1.sh {params.relatedness_file} {params.relatedness_thresh} {params.relatedness_file_filter};"

        "vcftools --vcf {input.filter1} --out {params.output_name} --remove {params.relatedness_file_filter} --recode --recode-INFO-all;"

#6) Correct for Missingness (SNP and Sample)
#I Perform SNP and Sample Missingness at the same time so that the process doesn't go on ad infinitum. Otherwise after each successful filter, I would have to perform the other filter. Read some papers to see how others implement this

rule Filter3:
    input:
        filter2 = interFolder + '/{sample}' + '_filter2.recode.vcf'
    output:
        filter3 = outFolder + '/{sample}' + '_QCFinished.recode.vcf'
    params:
        missingness_file_name = interFolder + '/{sample}',
        lmiss = interFolder + '/{sample}.lmiss',
        lmiss_filtered = interFolder + '/{sample}_lmiss_filtered.txt',
        imiss = interFolder + '/{sample}.imiss',
        imiss_filtered = interFolder + '/{sample}_imiss_filtered.txt',
        output_name = outFolder + '/{sample}' + '_QCFinished',
        miss_thresh_snp = MISS_THRESH_SNP,
        miss_thresh_indi = MISS_THRESH_INDI
    shell:
        "module load vcftools/0.1.15;"

        #Create Missingness files
        "vcftools --vcf {input.filter2} --out {params.missingness_file_name} --missing-indv;"
        "vcftools --vcf {input.filter2} --out {params.missingness_file_name} --missing-site;"

        #Create file of SNPs to remove from vcf due to missingness
        "./snakeRefScript1.sh {params.lmiss} {miss_thresh_snp} {params.lmiss_filtered};"

        #Create file of Individuals to remove from vcf due to missingness
        "./snakeRefScript2.sh {params.imiss} {miss_thresh_indi} {params.imiss_filtered};"

        #remove samples that do not pass threshold
        "vcftools --vcf {input.filter2} --out {params.output_name} --exclude-positions {params.lmiss_filtered} --remove {params.imiss_filtered} --recode --recode-INFO-all;"
