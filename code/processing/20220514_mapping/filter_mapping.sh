# Find working directiory
result=${PWD##*/}

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

# Go back path
parent_path=${parent_path%/*}
parent_path=${parent_path%/*}
parent_path=${parent_path%/*}

# Find data directory
folder=$parent_path'/data/barcodes/'$result'/per_gene/'



for FILE in $folder'*'
do
  awk -F"," '$3>2 {print $0}' $FILE > ${FILE%.*}'_filtered.txt'
done
