#!/bin/bash

# index for files
PARENT_PATH=$(dirname $(greadlink -f $0))
result=${PARENT_PATH##*/}

# Go back path
PARENT_PATH=${PARENT_PATH%/*}
PARENT_PATH=${PARENT_PATH%/*}
PARENT_PATH=${PARENT_PATH%/*}

# Find data directory
FOLDER=$PARENT_PATH'/data/filtered_sequencing/'$result
OUT_DIR=$PARENT_PATH'/data/extracted_barcodes/'$result
# Make directories for stored data
mkdir $OUT_DIR


# Sequence file names
file_bc=$FOLDER'/R1.fastq.gz'

# Find barcodes
gunzip -c $file_bc | awk ' NR%4==2 {
        print $0;
    }'| sort | uniq -c | sort -bgr > $OUT_DIR'/collapsed.txt'
