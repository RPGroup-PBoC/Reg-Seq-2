#!/bin/bash
# Group number
group=${1:-efflux}

# Find working directiory

PARENT_PATH=$(dirname $(greadlink -f $0))
result="${PARENT_PATH##*/}"
echo $result

# Go back path
PARENT_PATH=${PARENT_PATH%/*}
PARENT_PATH=${PARENT_PATH%/*}
PARENT_PATH=${PARENT_PATH%/*}

# Find data directory
processing_folder=$PARENT_PATH'/data/processed_promoters/'$result
sam_file=$processing_folder'/'$group'_collapsed.sam'

BBMAP_PATH=$(find $PARENT_PATH -name "bbmap.sh")

# Make directories for stored data
mkdir $PARENT_PATH'/data/barcodes/'$result
out_dir=$PARENT_PATH'/data/barcodes/'$result'/'$group'_per_gene'
if [ -d $out_dir ] 
then
    echo $out_dir" exists"
else 
    mkdir $out_dir
fi

echo $BBMAP_PATH

$BBMAP_PATH ref=$PARENT_PATH'/data/twist_orders/efflux_pumps_wt_sequences.fasta' path=$processing_folder

$BBMAP_PATH ambiguous='best' indelfilter='0' nfilter='0' minid='0.85' trimreaddescriptions='t' in=$processing_folder/$group'_collapsed.fasta' out=$sam_file t='8' path=$processing_folder


# Start the process
start=$SECONDS

# Input file name
echo $out_dir
echo "Assigning gene names to promoter-barcode pairs..."

# Extract promoter-bc pairs and corresponding gene names
awk -v o="$out_dir" 'BEGIN{FS="\t";OFS = ","} !(NR%500000){print NR " Promoters Processed"}; NF>10{gsub(/_/, ",", $1); print $10,$1,$3 >> (o"/"$3"_barcodes.txt")}' "$sam_file"

# terminal output message
echo "All Promoters Processed, now adding headers..."

# Loop through output directory and add a header to each file
cd $out_dir
echo "promoter,barcode,counts,name" > headerfile
for file in *barcodes.txt; do cat headerfile $file > tmpfile2; mv tmpfile2 "$file"; done

rm headerfile

# terminal output message
echo "done! Output files written to " "$out_dir"
end=$SECONDS
duration=$(( end - start ))
echo
echo "time elapsed: $duration seconds"
