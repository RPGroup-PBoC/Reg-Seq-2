N=${1:-1}
result=${PWD##*/}
parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

parent_path=${parent_path%/*}
parent_path=${parent_path%/*}
parent_path=${parent_path%/*}
folder=$parent_path'/data/processed_sequencing/'$result
mkdir $parent_path'/data/processed_promoters/'
mkdir $parent_path'/data/processed_promoters/'$result
cd $folder
FILES=$(find . -name '*R2*.fastq.gz')
(
for file in $FILES
do
  ((i=i%N)); ((i++==0)) && wait
  FILENAME=${file%.fastq.gz}
  FILENAME=${FILENAME#*/}
  echo $FILENAME
  gunzip -c $file |  awk 'NR%4==2 {print substr($0,0,160)}' > "$parent_path/data/processed_promoters/$result/${FILENAME}_promoters.txt" &
done
)
