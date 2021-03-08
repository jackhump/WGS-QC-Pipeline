# script for downloading and assembling the Ensembl dbSNP reference VCF

ml axel
ml bcftools

# download each file
for i in {1..22}; do 
    if [ ! -f homo_sapiens-chr${i}.vcf.gz ]; then
    	axel ftp://ftp.ensembl.org/pub/release-99/variation/vcf/homo_sapiens/homo_sapiens-chr${i}.vcf.gz
    	tabix homo_sapiens-chr${i}.vcf.gz
    fi
done

# set chr names
for i in {1..22}; do 
	echo $i chr$i 
done > chr_name_conv.txt

# merge together and set chr names to chr1
bcftools concat --threads 4 homo_sapiens-chr{1..22}.vcf.gz -Oz -o ensembl_v99_hg38_nochr.vcf.gz

bcftools annotate --rename-chrs chr_name_conv.txt -Oz -o ensembl_v99_hg38.vcf.gz ensembl_v99_hg38_nochr.vcf.gz

tabix -p vcf ensembl_v99_hg38.vcf.gz

rm -rf homo_sapiens-chr*