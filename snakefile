#0) Introductions

#snakemake pipeline for QC of vcf files. Created first on 191115.

#Debugging
#Error1: config.yaml not referenced properly. Seems that if statement to input default config.yaml is not accepting the default file. Come back
#Error2: make sure to delete Cluster folder on subsequent runs or else it will cause an error. Rename Cluster to include runID
#Error3: config.yaml format is wrong. No "=", instead ": "
#Error4: All snakemake shell scripts must be in quotes. Create a separate script to reference commands so I don't have to escape quotes every time
#Error5: None of the inputs, outputs, or params should have $ in front of variables
#Error6: Commas in lists of inputs, outputs, or params
#Error7: gzvcf_suffix not defined. Added in the line to define it.
#Error8: Need to create a rule all to define which files to output
#Error9: Correct rule all so that it uses expand function
#Error10: vcftools --out has file extension. Will declare parameter with filenames to make output file
#Error11: Error in rule1 code. Tested it in command line and it works fine. Best idea is that module loaded before is not retained
#Error12: Error11's fix might be necessary, but also, need to put semi colons at the end of each bash line or it doesnt know to end apparently. Idk if this expands to the scripts run as well, but i dont think so
#Error13: Missed a "params." in 2 places in rule 2
#Error14: It seems like snakemake does not like how i executed a script within the shell section. Try ./ instead of source. I don't want to have to escape every quote
#Error15: Jk Error14 is from snakeRefScript1 itself. Also I accidentally call the wrong script in rule Filter3. Fixed the second point.

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
        output_name = interFolder + '/{sample}' + '_filter2'
    shell:
        "module load vcftools/0.1.15;"
        "vcftools --vcf {input.filter1} --out {params.relatedness_file_name} --relatedness2;"

        "./snakeRefScript1.sh {params.relatedness_file} $RELATEDNESS_THRESH {params.relatedness_file_filter};"

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
        output_name = outFolder + '/{sample}' + '_QCFinished'
    shell:
        "module load vcftools/0.1.15;"

        #Create Missingness files
        "vcftools --vcf {input.filter2} --out {params.missingness_file_name} --missing-indv;"
        "vcftools --vcf {input.filter2} --out {params.missingness_file_name} --missing-site;"

        #Create file of SNPs to remove from vcf due to missingness
        #Double check the output file versus the actual missingness in file just in case its picking up the wrong lines
        "./snakeRefScript1.sh {params.lmiss} $MISS_THRESH_SNP {params.lmiss_filtered};"

        #Create file of Individuals to remove from vcf due to missingness
        "./snakeRefScript2.sh {params.imiss} $MISS_THRESH_INDI {params.imiss_filtered};"

        #remove samples that do not pass threshold
        "vcftools --vcf {input.filter2} --out {params.output_name} --exclude-positions {params.lmiss_filtered} --remove {params.imiss_filtered} --recode --recode-INFO-all;"
