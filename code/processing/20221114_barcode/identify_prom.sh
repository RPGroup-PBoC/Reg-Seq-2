gc=$1

PARENT_PATH=$(dirname $(greadlink -f $0))
RESULT=${PARENT_PATH##*/}

# Go back path
PARENT_PATH=${PARENT_PATH%/*}
PARENT_PATH=${PARENT_PATH%/*}
PARENT_PATH=${PARENT_PATH%/*}

# Data FOLDER
MAPPED=$PARENT_PATH'/data/barcodes/20220514_mapping/mapped_barcodes.csv'
FOLDER=$PARENT_PATH'/data/processed_barcodes/'$RESULT


awk 'NR==FNR{a[$1]=$0;next} ($3 in a){b=$3;$3="";print $0 a[b]}'  $MAPPED $FOLDER"/$(gc)_combined.txt" > $FOLDER"/$(gc)_identified.txt"
echo -e "cDNA_count gDNA_count barcode name mapping_count promoter" | cat - $FOLDER"/$(gc)_identified.txt" > "$FOLDER/$(gc)_identified_.txt"