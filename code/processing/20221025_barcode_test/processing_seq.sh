#!/bin/bash
# Find working directiory
PARENT_PATH=$(dirname $(greadlink -f $0))
RESULT=${PARENT_PATH##*/}

# Go back path
PARENT_PATH=${PARENT_PATH%/*}
PARENT_PATH=${PARENT_PATH%/*}
PARENT_PATH=${PARENT_PATH%/*}

# Data FOLDER
FOLDER=$PARENT_PATH'/data/sequencing/'$RESULT
DATA_FOLDER=$FOLDER"/*"
OUT_FOLDER=$PARENT_PATH'/data/processed_sequencing/'$RESULT

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
  OUT1=$OUT_FOLDER"/"$GC"_"$XNA"_R1.fastq.gz"
  HTML=$OUTPUT_DIR"/"$RESULT"/"$GC"_"$XNA"_fastp_report.html"
  JSON=$OUTPUT_DIR"/"$RESULT"/"$GC"_"$XNA"_fastp_report.json"
  fastp --in1 $FILE --out1 $OUT1 --trim_tail1 '6'  --verbose --disable_length_filtering --html $HTML --json $JSON --report_title $html_report --thread '6'
done
