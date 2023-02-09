#!/bin/bash
# Find working directiory
DATA_FOLDER="/home/tom/git/1000_genes_ecoli/data"
SRA=$1
OUT=$2
IN_PATH=$DATA_FOLDER"/sequencing/$SRA.fastq"
OUT_PATH=$DATA_FOLDER"/filtered_sequencing/"$OUT"_barcodes.fastq"

HTML=$DATA_FOLDER"/mapping_fastp_report.html"
JSON=$DATA_FOLDER"/mapping_fastp_report.json"
fastp -i $IN_PATH -o $OUT_PATH  --verbose --html $HTML --json $JSON --report_title $html_report --thread '12'  --average_qual '30'  --n_base_limit '0' --unqualified_percent_limit '10'
