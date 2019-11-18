#1) Input Variables
missingness=$1
cutoff=$2
output_file=$3

#2) Filter missingness based on cutoff and output in output_file
for i in $(seq 2 $(wc -l < $missingness))
do
if (( $(echo "$(sed -n "${i},${i}p" $missingness | cut -f 6) > $cutoff" | bc -l) ))
then echo -e "$(sed -n "${i},${i}p" $missingness | cut -f 1)" >> $output_file
fi
done

#Old Code
#for i in $(seq 2 $(wc -l < {params.imiss}))
#    do
#    if (( $(echo "$(sed -n "${i},${i}p" {params.imiss} | cut -f 6) > $MISS_THRESH_INDI" | bc -l) ))
#    then echo -e "$(sed -n "${i},${i}p" {params.imiss} | cut -f 1)" >> {params.imiss_filtered}
#    fi
#    done
