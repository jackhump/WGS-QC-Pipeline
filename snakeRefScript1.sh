#1) Input Variables
relatedness=$1
cutoff=$2
output_file=$3

#2) Filter relatedness based on cutoff and output in output_file
for i in $(seq 2 $(wc -l < $relatedness))
do
if (( $(echo "$(sed -n "${i},${i}p" $relatedness | cut -f 7) > $cutoff" | bc -l) ))
then if [[ "$(sed -n "${i},${i}p" $relatedness | cut -f 1)" != "$(sed -n "${i},${i}p" $relatedness | cut -f 2)" ]]
then echo -e "$(sed -n "${i},${i}p" $relatedness | cut -f 1)" >> $output_file
fi
fi
done


#Old Code
    #for i in $(seq 2 $(wc -l < {params.relatedness_file}))
    #do
    #if (( $(echo "$(sed -n "${i},${i}p" {params.relatedness_file} | cut -f 7) > $RELATEDNESS_THRESH" | bc -l) ))
    #if [[ "$(sed -n "${i},${i}p" {params.relatedness_file} | cut -f 1)" != "$(sed -n "${i},${i}p" {params.relatedness_file} | cut -f 2)" ]]
    #then echo -e "$(sed -n "${i},${i}p" {params.relatedness_file} | cut -f 1)" >> {params.relatedness_file_filter}
    #fi
    #fi
    #done
