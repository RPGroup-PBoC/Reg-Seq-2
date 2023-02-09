#!/bin/bash

# File containing barcode-promoter key
IN_MAP="/Volumes/rp_lab_ext/og_regseq_data/mapping_identified.csv"
# File containing barcode counts
IN_DNA_BC="/Volumes/rp_lab_ext/og_regseq_data/bc_by_gc/"${1:-1}"_DNA.txt"
IN_RNA_BC="/Volumes/rp_lab_ext/og_regseq_data/bc_by_gc/"${1:-1}"_RNA.txt"
# location to store identified barcodes
OUT_BC="/Volumes/rp_lab_ext/og_regseq_data/bc_by_gc/"${1:-1}"_identified.txt"

echo $IN_MAP
echo $IN_DNA_BC
awk 'NR==FNR{a[$1]=$0; next} ($1 in a){b=$1;$1="";print $0" "a[b]}' $IN_MAP $IN_DNA_BC > "temp_DNA.txt"
awk 'NR==FNR{a[$1]=$0; next} ($1 in a){b=$1;$1="";print $0" "a[b]}' $IN_MAP $IN_RNA_BC > "temp_RNA.txt"

awk 'NR==FNR{a[$3]=$0; next} ($3 in a){b=$3;print $2" "a[b]}' "temp_DNA.txt" "temp_RNA.txt" > "temp.txt"

echo -e "ct_1 index ct_0 barcode promoter mapping_count name nmut "| cat - "temp.txt"  > $OUT_BC

rm "temp.txt"
rm "temp_DNA.txt"
rm "temp_RNA.txt"