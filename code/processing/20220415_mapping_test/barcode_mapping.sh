N=${1:-1}
result=${PWD##*/}
parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

parent_path=${parent_path%/*}
parent_path=${parent_path%/*}
parent_path=${parent_path%/*}
folder=$parent_path'/data/processed_sequencing/'$result
mkdir $parent_path'/data/barcodes/'$result
cd $folder
FILES=$(find . -name '*merged.fastq.gz')
(
for file in $FILES
do
  ((i=i%N)); ((i++==0)) && wait
  FILENAME=${file%.fastq.gz}
  FILENAME=${FILENAME#*/}
  echo $FILENAME
  gunzip -c $file | awk /CGGTTTATGGGTGTTATCGC/ | awk '{ if ((length($0) == 272) || (length($0) == 252) || (length($0) == 233)) { bc = substr($0,0,20); prom=substr($0,length($0)-159,160); print bc "," prom } }' > "$parent_path/data/barcodes/$result/${FILENAME}_barcodes.txt" &
done
)
