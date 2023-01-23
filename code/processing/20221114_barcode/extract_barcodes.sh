#!/bin/bash
PARENT_PATH=$(dirname $(greadlink -f $0))
RESULT=${PARENT_PATH##*/}

# Go back path
PARENT_PATH=${PARENT_PATH%/*}
PARENT_PATH=${PARENT_PATH%/*}
PARENT_PATH=${PARENT_PATH%/*}

# Find data directory
FOLDER=$PARENT_PATH'/data/filtered_sequencing/'$RESULT
DATA_FOLDER=$FOLDER"/*"
# Make directories for stored data
OUT_FOLDER=$PARENT_PATH'/data/extracted_barcodes/'$RESULT

echo $DATA_FOLDER

if [ -d $OUT_FOLDER ] 
then
    echo $OUT_FOLDER" exists"
else 
    mkdir $OUT_FOLDER
fi

for FILE in $DATA_FOLDER;do
  filename="${FILE##*/}"
  filename="${filename%.*}"
  tmp=${filename#*_}
  XNA=${tmp%%_*}
  GC="${filename%%_*}"
  # Find barcodes
  gunzip -c $FILE | awk ' NR%4==2 {
        print $0;
    }
    NR%4==0 {
        print $0;
    }
    ' > $OUT_FOLDER"/"$GC'_'$XNA'_barcodes.txt'

  # Sort and count unique combinations
  cat $OUT_FOLDER"/"$GC'_'$XNA'_barcodes.txt' | awk 'NR%2==1 {print $0}' | sort | uniq -c | sort -bgr > $OUT_FOLDER"/"$GC'_'$XNA'_collapsed.txt'
done