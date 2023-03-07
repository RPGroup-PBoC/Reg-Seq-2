#!/bin/bash

# Find working directiory
PARENT_PATH=$(dirname $(greadlink -f $0))
RESULT=${PARENT_PATH##*/}

# Go back path
PARENT_PATH=${PARENT_PATH%/*}
PARENT_PATH=${PARENT_PATH%/*}
PARENT_PATH=${PARENT_PATH%/*}

# Data FOLDER
FOLDER=$PARENT_PATH'/data/processed_sequencing/'$RESULT
OUT_FOLDER=$PARENT_PATH'/data/barcodes/'$RESULT

if [ -d $OUT_FOLDER ] 
then
    echo $OUT_FOLDER" exists"
else 
    mkdir $OUT_FOLDER
fi

FILE='scramble_merged.fastq.gz'

FILENAME=${FILE%.fastq.gz}
FILENAME=${FILENAME#*/}
echo $FILENAME
gunzip -c $FOLDER"/"$FILE | awk /CACGCCAGTTGTGAACATAA/ | awk '{ if (length($0) == 222) { bc = substr($0,length($0)-19,20); prom=substr($0,0,160); print bc "," prom } }' | sort | uniq -c | awk -v OFS="\t" '$1=$1' |  sort -bgr > "${OUT_FOLDER}/${FILENAME}_barcodes_.txt" 

echo -e "counts\tbarcode\tpromoter" | cat - "${OUT_FOLDER}/${FILENAME}_barcodes_.txt"  > "${OUT_FOLDER}/${FILENAME}_barcodes.txt" 
rm "${OUT_FOLDER}/${FILENAME}_barcodes_.txt" 