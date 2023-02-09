#!/bin/bash

# index for files
group=${1:-efflux}
# Find working directiory


PARENT_PATH=$(dirname $(greadlink -f $0))

# Go back path
PARENT_PATH=${PARENT_PATH%/*}
PARENT_PATH=${PARENT_PATH%/*}
PARENT_PATH=${PARENT_PATH%/*}

# Find data directory
file=$PARENT_PATH'/data/twist_orders/TWIST_sequences_niko_30000.csv'

cat $file | awk -F',' 'NR%1501==1 {print ">"substr($1, 2, 4); print substr($2, 21, 160)"\n";}' > $PARENT_PATH'/data/twist_orders/efflux_pumps_wt_sequences.fasta'

