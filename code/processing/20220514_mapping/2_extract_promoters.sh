#!/bin/bash

# index for files
group=${1:-110}


PARENT_PATH=$(dirname $(greadlink -f $0))
result=${PARENT_PATH##*/}

# Go back path
PARENT_PATH=${PARENT_PATH%/*}
PARENT_PATH=${PARENT_PATH%/*}
PARENT_PATH=${PARENT_PATH%/*}

#PARENT_PATH="/Volumes/rp_lab_ext/1000_genes_ecoli"

# Find data directory
INFOLDER=$PARENT_PATH'/data/filtered_sequencing/'$result/

# Make directories for stored data
mkdir $PARENT_PATH'/data/extracted_pairs/'
mkdir $PARENT_PATH'/data/extracted_pairs/'$result
OUT_FOLDER=$PARENT_PATH'/data/extracted_pairs/'$result


rm -rf $OUT_FOLDER'/temp/'
mkdir $OUT_FOLDER'/temp/'

command=paste
for i in $INFOLDER*$group*.gz; do
    command="$command <(gunzip -cd $i)"
done
echo $command

# split 
eval $command | awk 'NR%4==2 {print $0}' > $OUT_FOLDER'/temp_file.txt'
split -n 40 $OUT_FOLDER'/temp_file.txt' $OUT_FOLDER'/temp/small_chunk'
rm $OUT_FOLDER'/temp_file.txt'
for X in $OUT_FOLDER/temp/small_chunk*; do 
  filename="${X##*/}"
  sort < $X > $OUT_FOLDER'/temp/sorted-'$filename;
done

sort -m $OUT_FOLDER/temp/sorted-small_chunk* | uniq -c | awk '{
    print ">"$3"_"$1"\n"$2;
    print "";
}' > $OUT_FOLDER'/'$group'_collapsed.fasta'


rm -rf $OUT_FOLDER'/temp/'

#eval $command | awk 'NR%4==2 {print $0}' | sort --parallel='4' | uniq -c | sort --parallel='4' -bgr | awk '{
 #   print $3"\t"$2;
#}' > $group'_collapsed.txt'