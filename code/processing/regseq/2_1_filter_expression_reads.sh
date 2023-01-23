#!/bin/bash
# Find working directiory
DATA_FOLDER="/Volumes/rp_lab_ext/og_regseq_data"
SRA=$1
OUT=$2
IN_PATH="/Volumes/rp_lab_ext/og_regseq_data/raw_sequencing/$SRA.fastq"
OUT_PATH="/Volumes/rp_lab_ext/og_regseq_data/filtered_sequencing/"$OUT"_barcodes.fastq"

HTML=$DATA_FOLDER"/mapping_fastp_report.html"
JSON=$DATA_FOLDER"/mapping_fastp_report.json"
fastp -i $IN_PATH -o $OUT_PATH  --verbose --html $HTML --json $JSON --report_title $html_report --thread '6'  --average_qual '33'  --n_base_limit '0' --unqualified_percent_limit '10'
