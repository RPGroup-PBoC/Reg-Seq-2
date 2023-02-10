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
DATA_FOLDER=$FOLDER"/230209_M05340_0411_000000000-KP6G8/Data/Intensities/BaseCalls"
OUT_FOLDER=$PARENT_PATH'/data/filtered_sequencing/'$RESULT

if [ -d $OUT_FOLDER ] 
then
    echo $OUT_FOLDER" exists"
else 
    mkdir $OUT_FOLDER
fi

FILES=("RP88" "RG88" "RP260" "RG260" "DG" "DP")

for FILE in "${FILES[@]}";do
for _FILE in $DATA_FOLDER"/"$FILE*R1*;do
  echo $_FILE
  OUT1=$OUT_FOLDER"/"$FILE".fastq.gz"
  HTML=$OUT_FOLDER"/"$FILE"_fastp_report.html"
  JSON=$OUT_FOLDER"/"$FILE"_fastp_report.json"
  fastp --in1 $_FILE --out1 $OUT1 --trim_tail1 '6'  --verbose --disable_length_filtering --html $HTML --json $JSON --report_title $html_report --thread '6' -q '20' --n_base_limit '0' --unqualified_percent_limit '10'
done
done
