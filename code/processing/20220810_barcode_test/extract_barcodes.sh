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
folder=$PARENT_PATH'/data/processed_sequencing/'$result

# Make directories for stored data
mkdir $PARENT_PATH'/data/processed_barcodes/'
mkdir $PARENT_PATH'/data/processed_barcodes/'$result

# Go to data
cd $folder

# Sequence file names
file_bc=$group'_R1.fastq.gz'

# Find barcodes
gunzip -c $file_bc | awk ' NR%4==2 {
        print $0;
    }
    NR%4==0 {
        print $0;
    }
    ' > $group'_barcodes.txt'


# Sort and count unique combinations
cat $group'_barcodes.txt' | awk 'NR%2==1 {print $0}' | sort | uniq -c | sort -bgr > $group'_collapsed.txt'


mv $group'_collapsed.txt' $PARENT_PATH'/data/processed_barcodes/'$result/$group'_collapsed.txt'