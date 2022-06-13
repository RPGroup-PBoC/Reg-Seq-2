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
folder=$parent_path'/data/processed_promoters/'$result

# make key
../../../bbmap/bbmap.sh ref=$parent_path'/data/wt_sequences.fasta'

../../../bbmap/bbmap.sh ambiguous='best' indelfilter='0' nfilter='0' minid='0.85' trimreaddescriptions='t' in=$folder/$group'_collapsed.fasta' out=$folder/$group'_collapsed.sam' t=8