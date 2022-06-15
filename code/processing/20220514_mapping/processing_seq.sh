#!/bin/bash
GROUP=${1:-110}
# Find working directiory
RESULT=${PWD##*/}

PARENT_PATH=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

# Go back path
PARENT_PATH=${PARENT_PATH%/*}
PARENT_PATH=${PARENT_PATH%/*}
PARENT_PATH=${PARENT_PATH%/*}

# Data FOLDER
FOLDER=$PARENT_PATH'/data/sequencing/'$RESULT
OUT_FOLDER=$PARENT_PATH'/data/processed_sequencing/'$RESULT

if [ -d $OUT_FOLDER ] 
then
    echo $OUT_FOLDER" exists"
else 
    mkdir $OUT_FOLDER
fi


# Find read1 and read2 files
READ1=$(find $FOLDER -name $GROUP"*R1*.fastq.gz")
READ2=$(find $FOLDER -name $GROUP"*R2*.fastq.gz")

# Define output file paths
OUT1=$OUT_FOLDER"/"$GROUP"_R1.fastq.gz"
OUT2=$OUT_FOLDER"/"$GROUP"_R2.fastq.gz"

HTML=$OUTPUT_DIR"/"$RESULT"/"$GROUP"_fastp_report.html"
JSON=$OUTPUT_DIR"/"$RESULT"/"$GROUP"_fastp_report.json"

# Define string to be ran on the termina
fastp --in1 $READ1 --in2 $READ2 --out1 $OUT1 --out2 $OUT2 --trim_tail1 '11' --trim_tail2 '11' --verbose --disable_length_filtering --html $HTML --json $JSON --report_title $html_report --thread '6'