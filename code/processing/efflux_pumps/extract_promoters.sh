#!/bin/bash

# index for files
group=${1:-efflux}
# Find working directiory


PARENT_PATH=$(dirname $(greadlink -f $0))
result=${PARENT_PATH##*/}

# Go back path
PARENT_PATH=${PARENT_PATH%/*}
PARENT_PATH=${PARENT_PATH%/*}
PARENT_PATH=${PARENT_PATH%/*}

# Find data directory
folder=$PARENT_PATH'/data/processed_sequencing/'$result

# Make directories for stored data
mkdir $PARENT_PATH'/data/processed_promoters/'
mkdir $PARENT_PATH'/data/processed_promoters/'$result

# Go to data
cd $folder

# Sequence file names
file=$group'_merged_rev.fastq.gz'
# Find barcodes
gunzip -c $file | awk ' NR%4==2 {
        print substr($1, 21, 160)"\t"substr($1, 221, 20);
    }
    ' > $group'_combined.txt'

# Sort and count unique combinations
cat $group'_combined.txt' | sort | uniq -c | sort -bgr > $group'_collapsed.txt'

# Remove temporary files
rm $group'_combined.txt'

awk '{
    print ">"$3"_"$1"\n"$2;
    print "";
}'  $group'_collapsed.txt' > $group'_collapsed.fasta'

mv $group'_collapsed.fasta' $PARENT_PATH'/data/processed_promoters/'$result/$group'_collapsed.fasta'

TTTTCTACTTTCCGGCTTGCATAGTTATGAATACAGGCATCTCAAGGCACATAAACACAAAAAAAGATTAATATTCTACTGTTTTATTTTGACGCGGGTTGAAAGAGGCAGAATTAAAACCTCGTAAATTGAAATATATATTGATGTAGTGAATGTATCTTAGGTAAATAATATATATTAATGGTAAGAAGCTCCCACAA