#!/bin/bash
PARENT_PATH=$(dirname $(greadlink -f $0))
RESULT=${PARENT_PATH##*/}

# Go back path
PARENT_PATH=${PARENT_PATH%/*}
PARENT_PATH=${PARENT_PATH%/*}
PARENT_PATH=${PARENT_PATH%/*}
PARENT_PATH="/Volumes/rp_lab_ext/1000_genes_ecoli"

# Find data directory
FOLDER=$PARENT_PATH'/data/filtered_sequencing/'$RESULT
DATA_FOLDER=$FOLDER"/*"
# Make directories for stored data
OUT_FOLDER=$PARENT_PATH'/data/extracted_barcodes/'$RESULT


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
  echo "processing "$XNA" for "$GC
  echo
  # Find barcodes
  gunzip -c $FILE | awk ' NR%4==2 {
        print $0;
    }
    ' > $OUT_FOLDER"/"$GC'_'$XNA'_barcodes.txt'
  echo "Barcodes extracted..."
  echo
  # Sort and count unique combinations
  echo "Splitting file"
  echo
  split -n 40 $OUT_FOLDER"/"$GC'_'$XNA'_barcodes.txt' $OUT_FOLDER'/temp/small_chunk' 
  echo "Sorting individual files"
  echo
  for X in $OUT_FOLDER/temp/small_chunk*; do 
    filename="${X##*/}"
    sort < $X > $OUT_FOLDER'/temp/sorted-'$filename;
  done
  echo "Merging files and counting unique barcodes"
  echo
  sort -m $OUT_FOLDER/temp/sorted-small_chunk* | uniq -c > $OUT_FOLDER"/"$GC'_'$XNA'_collapsed.txt'
  rm $OUT_FOLDER/temp/sorted-*
  rm $OUT_FOLDER/temp/small_chunk*
  echo "DONE"
  echo
done