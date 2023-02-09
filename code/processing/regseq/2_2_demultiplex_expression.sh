#!/bin/bash

IN=$1

<<<<<<< HEAD
DATA_FOLDER="/home/tom/git/1000_genes_ecoli/data"
IN_PATH="/home/tom/git/1000_genes_ecoli/data/filtered_sequencing/"$IN"_barcodes.fastq"
rm -rf $DATA_FOLDER"/"$IN"_bc_by_gc/"
mkdir $DATA_FOLDER"/"$IN"_bc_by_gc/"

rm -rf $DATA_FOLDER"/"$IN"_bc_by_gc/temp/"
mkdir $DATA_FOLDER"/"$IN"_bc_by_gc/temp/"


declare -A dict_group

while IFS=' ' read -r value key; do
    dict_group[$key]=$value
done < $DATA_FOLDER/index_group.txt


declare -A dict_gc

while IFS=' ' read -r key value; do
    dict_gc[$key]=$value
done < $DATA_FOLDER/gc_barcodes.txt

echo "Filtering..."
cat $IN_PATH | awk '/TATTAGGCTTCTCCTCAGCG/' | awk '/TCACTGGCCGTCGTTTTACATGACTGACTGA/' | awk -v o=$DATA_FOLDER"/"$IN"_bc_by_gc/temp" 'FNR==1{++f} \
f==1 {a[$2]=$1} \
f==2 {b[$1]=$2} \
f==3 {ind1=substr($0, 0, 4); bc=substr($0,60,20); group=substr($0, 25, 4); if((group in a) && (ind1 in b)){
printf "%s\n", bc >>o"/"b[ind1]"_"a[group]".txt"}}' $DATA_FOLDER/index_group.txt $DATA_FOLDER/gc_barcodes.txt -


echo "Counting unique barcodes..."
OUT=$DATA_FOLDER"/"$IN"_bc_by_gc"/temp/*.txt
for FILE in $OUT;do
    filename="${FILE##*/}"
    filename="${filename%.*}"
    sort --parallel 20 -T ./ $FILE | uniq -c | sort --parallel 20  -bgr -T ./|  awk -v OFS="\t" '$1=$1' > $DATA_FOLDER"/"$IN"_bc_by_gc/"$filename"_expression_counted.csv";
done
=======
DATA_FOLDER="/Volumes/rp_lab_ext/og_regseq_data"
IN_PATH="/Volumes/rp_lab_ext/og_regseq_data/filtered_sequencing/"$IN"_barcodes.fastq"
rm -rf $DATA_FOLDER"/"$IN"_bc_by_index/"
mkdir $DATA_FOLDER"/"$IN"_bc_by_index/"


cat $IN_PATH | awk '/TATTAGGCTTCTCCTCAGCG/' | awk '/TCACTGGCCGTCGTTTTACATGACTGACTGA/' | awk '{ind1=substr($0, 0, 4); ind2=substr($0, 25, 4); bc=substr($0,60,20); printf "%s\t%s\t%s\n", ind1, ind2, bc}' | sort --parallel 6 | uniq -c | sort -bgr |  awk -v OFS="\t" '$1=$1' > $DATA_FOLDER"/"$IN"_expression_counted.csv"

awk -v o="$DATA_FOLDER/"$IN"_bc_by_index" '{if(length($2)==4) {print $1,$3,$4 >> (o"/"$2".txt")}}' $DATA_FOLDER"/"$IN"_expression_counted.csv"
>>>>>>> 0994e8cb042ecb827a60dd9eb06d37c7a5029da8
