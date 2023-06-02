#!/bin/bash

# index for files
group=${1:-1}
x=$2
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
file_bc=$FOLDER/$group'_R1.fastq.gz'

# Compute number of lines and how much to sample
N=$(gunzip -c $file_bc | sed -n '$=')
N_INT=$(($N / $x / 4))
# Find barcodes
    
gunzip -c $file_bc | awk ' NR%4==2 {
        print $0;
    }'| shuf -n "${N_INT}" | sort | uniq -c | sort -bgr > $OUT_DIR/$group'_collapsed_shuff_'$x'.txt'
