#!/bin/bash

# Find working directiory
PARENT_PATH=$(dirname $(greadlink -f $0))
RESULT=${PARENT_PATH##*/}

# Go back path
PARENT_PATH=${PARENT_PATH%/*}
PARENT_PATH=${PARENT_PATH%/*}
PARENT_PATH=${PARENT_PATH%/*}


# Data FOLDER
FOLDER=$PARENT_PATH'/data/sequencing/'$RESULT
OUT_FOLDER=$PARENT_PATH'/data/filtered_sequencing/'$RESULT

if [ -d $OUT_FOLDER ] 
then
    echo $OUT_FOLDER" exists"
else 
    mkdir $OUT_FOLDER
fi


# Find read1 and read2 files
READ1=$(find $FOLDER -name "*Undetermined*R1*.fastq.gz")
READ2=$(find $FOLDER -name "*Undetermined*R2*.fastq.gz")

echo $READ1
echo $READ2

# Define output file paths
OUT1=$OUT_FOLDER"/R1.fastq.gz"
OUT2=$OUT_FOLDER"/R2.fastq.gz"

HTML=$OUT_FOLDER"/"$GROUP"_fastp_report.html"
JSON=$OUT_FOLDER"/"$GROUP"_fastp_report.json"

# Define string to be ran on the terminal
fastp --in1 $READ1 --in2 $READ2 --out1 $OUT1 --out2 $OUT2 --trim_tail1 '6' --trim_tail2 '6' --verbose --disable_length_filtering --html $HTML --json $JSON --thread '6' -q '20' --n_base_limit '0' --unqualified_percent_limit '10'