#!/bin/bash

# Find working directiory
DATA_FOLDER="/Volumes/rp_lab_ext/og_regseq_data"
IN="/Volumes/rp_lab_ext/og_regseq_data/filtered_sequencing/mapping_merged.fastq"

cat $IN | awk 'NR%4==2 { if (length($0) == 295 || length($0)==299){ l=length($0)-295;prom=substr($0, 21, 160); bc=substr($0,256,20); printf "%s\t%s\n", prom, bc};}' | sort | uniq -c | sort -bgr |  awk -v OFS="\t" '$1=$1' > $DATA_FOLDER"/mapping_counted.csv"