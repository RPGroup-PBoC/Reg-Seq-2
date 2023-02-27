#!/bin/bash

# index for files

PARENT_PATH=$(dirname $(greadlink -f $0))
RESULT=${PARENT_PATH##*/}

# Go back path
PARENT_PATH=${PARENT_PATH%/*}
PARENT_PATH=${PARENT_PATH%/*}
PARENT_PATH=${PARENT_PATH%/*}

# 
DATA_FOLDER=$PARENT_PATH'/data/filtered_sequencing/'$RESULT/
OUT_FOLDER=$PARENT_PATH'/data/extracted_pairs/'$RESULT
# Make directories for stored data
mkdir $PARENT_PATH'/data/extracted_pairs/'
mkdir $PARENT_PATH'/data/extracted_pairs/'$RESULT


command=paste
for i in $DATA_FOLDER*.gz; do
    command="$command <(gunzip -cd $i)"
done
echo $command

# split 
eval $command | awk 'NR%4==2 {print $0}' | sort --parallel='4' | uniq -c | sort --parallel='4' -bgr | awk '{
 print ">"$3"_"$1"\n"$2;
    print "";
}' > $OUT_FOLDER/'collapsed.fasta'

