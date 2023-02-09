#!/bin/bash


gc=${1:-1}
FOLDER="/Volumes/rp_lab_ext/og_regseq_data/bc_by_gc/"
IN_DNA=$FOLDER$gc"_DNA_identified.txt"
IN_RNA=$FOLDER$gc"_RNA_identified.txt"


awk 'NR==FNR{a[$2]=$0; next} ($2 in a){print $1" "a[$2]}' $IN_DNA $IN_RNA | sort -bgr > "temp.txt"
echo -e "ct_1 ct_0 barcode name mapping_count promoter nmut " | cat - "temp.txt"  > $FOLDER$gc"_merged.txt"

rm "temp.txt"

