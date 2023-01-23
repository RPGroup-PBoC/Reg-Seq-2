#!/bin/bash

IN=$1

DATA_FOLDER="/Volumes/rp_lab_ext/og_regseq_data"
IN_PATH="/Volumes/rp_lab_ext/og_regseq_data/filtered_sequencing/"$IN"_barcodes.fastq"
rm -rf $DATA_FOLDER"/"$IN"_bc_by_index/"
mkdir $DATA_FOLDER"/"$IN"_bc_by_index/"


cat $IN_PATH | awk '/TATTAGGCTTCTCCTCAGCG/' | awk '/TCACTGGCCGTCGTTTTACATGACTGACTGA/' | awk '{ind1=substr($0, 0, 4); ind2=substr($0, 25, 4); bc=substr($0,60,20); printf "%s\t%s\t%s\n", ind1, ind2, bc}' | sort --parallel 6 | uniq -c | sort -bgr |  awk -v OFS="\t" '$1=$1' > $DATA_FOLDER"/"$IN"_expression_counted.csv"

awk -v o="$DATA_FOLDER/"$IN"_bc_by_index" '{if(length($2)==4) {print $1,$3,$4 >> (o"/"$2".txt")}}' $DATA_FOLDER"/"$IN"_expression_counted.csv"
