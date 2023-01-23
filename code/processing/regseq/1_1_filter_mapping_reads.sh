#!/bin/bash
# Find working directiory
DATA_FOLDER="/Volumes/rp_lab_ext/og_regseq_data"

IN1="/Volumes/rp_lab_ext/og_regseq_data/SRR10971993_1.fastq"
IN2="/Volumes/rp_lab_ext/og_regseq_data/SRR10971993_2.fastq"
OUT="/Volumes/rp_lab_ext/og_regseq_data/mapping_merged.fastq"

HTML=$DATA_FOLDER"/mapping_fastp_report.html"
JSON=$DATA_FOLDER"/mapping_fastp_report.json"
fastp --in1 $IN1 --in2 $IN2 --merged_out $OUT  --verbose --html $HTML --json $JSON --report_title $html_report --thread '6' --merge --overlap_len_require '1' --qualified_quality_phred '33' --n_base_limit '0' 
