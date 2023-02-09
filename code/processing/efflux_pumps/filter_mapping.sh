#!/bin/bash
group=${1:-efflux}
# Find working directiory


PARENT_PATH=$(dirname $(greadlink -f $0))
result="${PARENT_PATH##*/}"
# Go back path
PARENT_PATH=${PARENT_PATH%/*}
PARENT_PATH=${PARENT_PATH%/*}
PARENT_PATH=${PARENT_PATH%/*}



# Make directories if not existent
mkdir $PARENT_PATH'/data/barcodes/'$result'/'$group'_per_gene_filtered/'
mkdir $PARENT_PATH'/data/barcodes/'$result'/'$group'_per_gene/'

out_folder=$PARENT_PATH'/data/barcodes/'$result'/'$group'_per_gene_filtered/'
data_folder=$PARENT_PATH'/data/barcodes/'$result'/'$group'_per_gene/*'

for FILE in $data_folder;do
  GENE="${FILE##*/}"
  GENE="${GENE%.*}"
  awk -F"," '$3>2 {print $0}' $FILE > $out_folder$GENE'_filtered.txt'
done
