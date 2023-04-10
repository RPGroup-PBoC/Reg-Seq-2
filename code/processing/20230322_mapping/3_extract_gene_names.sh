#!/bin/bash

############################################################
# Process the input options.                               #
############################################################
# Get the options
#while getopts ":hi:o:bn:" option; do
#   case $option in
#      h) # display Help
#         Help
#         exit;;
#      i) # Input file or directory
#         in_file=$OPTARG;;
#      o) # Output directory
 #        out_dir=$OPTARG;;
#      b) # input file is BAM format
 #        Name=$OPTARG;;
 #     n) # Enter a name
#         Name=$OPTARG;;     
 #    \?) # Invalid option
 #        echo "Error: Invalid option"
#         exit;;
#   esac
#done


# Find working directiory

PARENT_PATH=$(dirname $(greadlink -f $0))
RESULT="${PARENT_PATH##*/}"
echo $RESULT

# Go back path
PARENT_PATH=${PARENT_PATH%/*}
PARENT_PATH=${PARENT_PATH%/*}
PARENT_PATH=${PARENT_PATH%/*}

# Find data directory
PROCESSING_FOLDER=$PARENT_PATH'/data/extracted_pairs/'$RESULT
SAM_FILE=$PROCESSING_FOLDER'/collapsed.sam'

BBMAP_PATH=$(find $PARENT_PATH -name "bbmap.sh")

# Make directories for stored data
out_dir=$PARENT_PATH'/data/barcodes/'$RESULT
if [ -d $out_dir ] 
then
    echo $out_dir" exists"
else 
    mkdir $out_dir
fi

echo $BBMAP_PATH

$BBMAP_PATH ref=$PARENT_PATH'/data/wt_sequences.fasta' path=$PROCESSING_FOLDER
echo $PROCESSING_FOLDER/'collapsed.fasta'
$BBMAP_PATH ambiguous='best' indelfilter='0' nfilter='0' minid='0.85' trimreaddescriptions='t' in=$PROCESSING_FOLDER/'collapsed.fasta' out=$SAM_FILE t='8' path=$PROCESSING_FOLDER


# Start the process
start=$SECONDS

# Input file name
echo $out_dir
echo "Assigning gene names to promoter-barcode pairs..."

# Extract promoter-bc pairs and corresponding gene names
awk -v o="$out_dir" 'BEGIN{FS="\t";OFS = ","}; NF>10{gsub(/_/, ",", $1); print $10,$1,$3,$5}' > $out_dir'/mapping_identified.csv' "$SAM_FILE"

# terminal output message
#echo "All Promoters Processed, now adding headers..."
#
# Loop through output directory and add a header to each file
#cd $out_dir
#echo "promoter,barcode,counts,name" > headerfile
#for file in *barcodes.txt; do cat headerfile $file > tmpfile2; mv tmpfile2 "$file"; done

#rm headerfile

# terminal output message
echo "done! Output files written to " "$out_dir"
end=$SECONDS
duration=$(( end - start ))
echo
echo "time elapsed: $duration seconds"
