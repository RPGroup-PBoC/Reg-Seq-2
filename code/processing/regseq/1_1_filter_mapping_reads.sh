#!/bin/bash
# Find working directiory
DATA_FOLDER="/home/tom/git/1000_genes_ecoli/data"

IN1=$DATA_FOLDER"/sequencing/SRR10971993_1.fastq"
IN2=$DATA_FOLDER"/sequencing/SRR10971993_2.fastq"
OUT=$DATA_FOLDER"/filtered_sequencing/mapping_merged.fastq"

HTML=$DATA_FOLDER"/mapping_fastp_report.html"
JSON=$DATA_FOLDER"/mapping_fastp_report.json"
fastp --in1 $IN1 --in2 $IN2 --merged_out $OUT  --verbose --html $HTML --json $JSON --report_title $html_report --thread '6' --merge --overlap_len_require '3' --n_base_limit '0' 
