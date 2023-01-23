#!/bin/bash
# Find working directiory
DATA_FOLDER="/Volumes/rp_lab_ext/og_regseq_data"
IN=$DATA_FOLDER"/expression_counted.csv"

# Extract promoter-bc pairs and corresponding gene names
awk -v o="$DATA_FOLDER/split_files" '{if(($2=="AGAG") || ($2=="CAAG") || ($2=="TCTA") || ($2=="ATGC")){print $1,$3 >> (o"/"$2".txt")}}' $IN