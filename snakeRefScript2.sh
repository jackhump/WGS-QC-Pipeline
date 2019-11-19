#1) Input Variables
missingness=$1
cutoff=$2
output_file=$3

#2) Filter missingness based on cutoff and output in output_file
for i in $(seq 2 $(wc -l < $missingness))
do
if (( $(echo "$(sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g' <<< $(sed -n "${i},${i}p" $missingness | cut -f 6)) > $cutoff" | bc -l) ))
#if (( $(echo "$(sed -n "${i},${i}p" $missingness | cut -f 6) > $cutoff" | bc -l) ))
then echo -e "$(sed -n "${i},${i}p" $missingness | cut -f 1)"'\t'"$(sed -n "${i},${i}p" $missingness | cut -f 2)" >> $output_file
fi
done


#Old code
#for i in $(seq 2 $(wc -l < {params.lmiss}))
#    do
#    if (( $(echo "$(sed -n "${i},${i}p" {params.lmiss} | cut -f 6) > $MISS_THRESH_SNP" | bc -l) ))
#    then echo -e "$(sed -n "${i},${i}p" {params.lmiss} | cut -f 1)"'\t'"$(sed -n "${i},${i}p" {params.lmiss} | cut -f 2)" >> {params.lmiss_filtered}
#    fi
#    done
