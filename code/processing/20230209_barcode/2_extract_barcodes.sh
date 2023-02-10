#!/bin/bash
# Find working directiory
PARENT_PATH=$(dirname $(greadlink -f $0))
RESULT=${PARENT_PATH##*/}

# Go back path
PARENT_PATH=${PARENT_PATH%/*}
PARENT_PATH=${PARENT_PATH%/*}
PARENT_PATH=${PARENT_PATH%/*}


# Data FOLDER
DATA_FOLDER=$PARENT_PATH'/data/filtered_sequencing/'$RESULT
OUT_FOLDER=$PARENT_PATH'/data/extracted_barcodes/'$RESULT

if [ -d $OUT_FOLDER ] 
then
    echo $OUT_FOLDER" exists"
else 
    mkdir $OUT_FOLDER
fi

FILES=("RP88" "RG88" "RP260" "RG260" "DG" "DP")

for FILE in "${FILES[@]}";do
  FILEPATH=$DATA_FOLDER/$FILE.fastq
  cat $FILEPATH | awk 'NR%4==2 {print $0}' | sort | uniq -c | sort -bgr > $OUT_FOLDER/${FILE}_collapsed.txt
done
