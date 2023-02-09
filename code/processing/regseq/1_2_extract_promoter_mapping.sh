#!/bin/bash

# Find working directiory
DATA_FOLDER="/home/tom/git/1000_genes_ecoli/data"
IN=$DATA_FOLDER"/filtered_sequencing/mapping_merged.fastq"
OUT=$DATA_FOLDER"/mapping_sequences"

rm -rf $OUT
mkdir $OUT

rm -rf $OUT/temp
mkdir $OUT/temp

declare -A dict

while IFS=' ' read -r value key; do
    dict[$key]=$value
done < $DATA_FOLDER/index_group.txt


# load file > find primer sequence > check for length > find index, barcode and promoter > store by index
echo "Filtering and extracting barcodes"
cat $IN | awk '/TATTAGGCTTCTCCTCAGCG/' | awk -v o=$OUT/temp '{ if (length($0) == 295 || length($0)==299)
    {ind_loc=index($0, "TATTAGGCTTCTCCTCAGCG"); Index=substr($0, ind_loc+20, 4);prom=substr($0, 21, 160); bc=substr($0,256,20); printf "%s\t%s\n", prom, bc >> (o"/"Index".txt") };}' 

echo "Count unique pairs"
for FILE in $OUT/temp/*;do
    filename="${FILE##*/}"
    filename="${filename%.*}" 
    if ([[ ${dict[$filename]} ]])
    then 
        sort -T "." --parallel 20 $FILE| uniq -c | sort -T "." -bgr  --parallel 20|  awk -v OFS="\t" '$1=$1' > $OUT"/"${dict[$filename]}"_mapping_counted.csv";
    fi
done