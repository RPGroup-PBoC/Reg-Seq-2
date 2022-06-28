#!/bin/bash

# index for files
group=${1:-110}


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

command=paste
for i in *$group*.gz; do
    command="$command <(gunzip -cd $i)"
done
eval $command | awk 'NR%4==2 {print $0}' | sort | uniq -c | sort -bgr | awk '{
    print ">"$3"_"$1"\n"$2;
    print "";
}' > $group'_collapsed.fasta'
