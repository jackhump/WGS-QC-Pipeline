
input=$1

output=$2


mkdir $output
#ml bcftools

bcftools query -f '%DP\n' $input > $output/$output.DP.txt

bcftools query -f '%VQSLOD\n' $input > $output/$output.VQSLOD.txt

bcftools query -f '%MQ\n' $input > $output/$output.MQ.txt

Rscript plot_quality_metrics.R $output/$output

