#!/bin/bash

# index for files
group=${1:-110}


PARENT_PATH=$(dirname $(greadlink -f $0))
result=${PARENT_PATH##*/}

# Go back path
PARENT_PATH=${PARENT_PATH%/*}
PARENT_PATH=${PARENT_PATH%/*}
PARENT_PATH=${PARENT_PATH%/*}

PARENT_PATH="/Volumes/rp_lab_ext/1000_genes_ecoli"

# Find data directory
folder=$PARENT_PATH'/data/filtered_sequencing/'$result

# Make directories for stored data
mkdir $PARENT_PATH'/data/extracted_pairs/'
mkdir $PARENT_PATH'/data/extracted_pairs/'$result

# Go to data
cd $folder

command=paste
for i in *$group*.gz; do
    command="$command <(gunzip -cd $i)"
done
echo $command

# split 
eval $command | awk 'NR%4==2 {print $0}' > $PARENT_PATH'/data/extracted_pairs/'$result'/temp_file.txt'
split -n 40 $PARENT_PATH'/data/extracted_pairs/'$result'/temp_file.txt' $PARENT_PATH'/data/extracted_pairs/'$result'/small_chunk'
for X in $PARENT_PATH/data/extracted_pairs/$result/small_chunk*; do 
  filename="${X##*/}"
  sort < $X > $OUT_FOLDER'/temp/sorted-'$filename;
done

sort -m $PARENT_PATH/sorted-small-chunk* | uniq -c | awk '{
    print ">"$3"_"$1"\n"$2;
    print "";
}' > $group'_collapsed.fasta'




#eval $command | awk 'NR%4==2 {print $0}' | sort --parallel='4' | uniq -c | sort --parallel='4' -bgr | awk '{
 #   print $3"\t"$2;
#}' > $group'_collapsed.txt'