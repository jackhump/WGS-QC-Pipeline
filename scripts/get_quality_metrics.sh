#ml bcftools
#ml R/3.6.2

# VCF
input=$1
# OUTDIR
output_dir=$2
# OUTPUT_NAME
output_name=$3
# Avg WGS Coverage
coverage=$4

mkdir $output_dir

bcftools query -f '%DP\n' $input > $output_dir/$output_name.DP.txt

bcftools query -f '%VQSLOD\n' $input > $output_dir/$output_name.VQSLOD.txt

bcftools query -f '%MQ\n' $input > $output_dir/$output_name.MQ.txt

bcftools query -l $input > $output_dir/$output_name.samples.txt

Rscript plot_quality_metrics.R $output_dir/$output_name $coverage
