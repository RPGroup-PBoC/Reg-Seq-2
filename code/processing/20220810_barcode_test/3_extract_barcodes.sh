#!/bin/bash

# index for files
group=${1:-1}
PARENT_PATH=$(dirname $(greadlink -f $0))
result=${PARENT_PATH##*/}

# Go back path
PARENT_PATH=${PARENT_PATH%/*}
PARENT_PATH=${PARENT_PATH%/*}
PARENT_PATH=${PARENT_PATH%/*}

# Find data directory
FOLDER=$PARENT_PATH'/data/filtered_sequencing/'$result
OUT_DIR=$PARENT_PATH'/data/extracted_barcodes/'
# Make directories for stored data
mkdir $PARENT_PATH'/data/filtered_sequencing/'
mkdir $PARENT_PATH'/data/filtered_sequencing/'$result


# Sequence file names
file_bc=$FOLDER/$group'_R1.fastq.gz'

# Find barcodes
gunzip -c $file_bc | awk ' NR%4==2 {
        print $0;
    }'| sort | uniq -c | sort -bgr > $OUT_DIR$group'_collapsed.txt'
