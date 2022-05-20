#!/bin/bash

# index for files
group=$1

# Find working directiory
result=${PWD##*/}

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

# Go back path
parent_path=${parent_path%/*}
parent_path=${parent_path%/*}
parent_path=${parent_path%/*}

# Find data directory
folder=$parent_path'/data/processed_sequencing/'$result

# Make directories for stored data
mkdir $parent_path'/data/processed_promoters/'
mkdir $parent_path'/data/processed_promoters/'$result

# Go to data
cd $folder

# Sequence file names
file_bc=$group'_R2.fastq.gz'
file_prom=$group'_R1.fastq.gz'

# Find barcodes
gunzip -c $file_bc | awk ' NR%4==2 {
        print $0;
    }
    ' > $group'_barcodes.txt'

# Find promoters
gunzip -c $file_prom | awk ' NR%4==2 {
        print $0;
    }
    ' > $group'_promoters.txt'

# Combine promoters and barcodes
paste $group'_barcodes.txt' $group'_promoters.txt' > $group'_combined.txt'

# Sort and count unique combinations
sort $group'_combined.txt' | uniq -c | sort -bgr > $group'_collapsed.txt'

# Do some formatting
tr -s "\t" " " < $group'_collapsed.txt' > $group'_collapsed.txt'

# Remove temporary files
rm $group'_promoters.txt'
rm $group'_barcodes.txt'
rm $group'_combined.txt'



