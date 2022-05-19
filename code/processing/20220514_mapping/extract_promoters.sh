group=$2

result=${PWD##*/}
parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

parent_path=${parent_path%/*}
parent_path=${parent_path%/*}
parent_path=${parent_path%/*}
folder=$parent_path'/data/processed_sequencing/'$result
mkdir $parent_path'/data/processed_promoters/'
mkdir $parent_path'/data/processed_promoters/'$result
cd $folder
file_bc=$(find . -name $group'*R2*.fastq.gz')
file_prom=$(find . -name $group'*R1*.fastq.gz')
(
for file in $FILES
do
  ((i=i%N)); ((i++==0)) && wait
  FILENAME_bc=${file_bc%.fastq.gz}
  FILENAME_bc=${FILENAME#*/}
  FILENAME_prom=${file_prom%.fastq.gz}
  FILENAME_prom=${FILENAME#*/}
  gunzip -c $file |  awk 'NR%4==2 {print substr($0,0,160)}' > "$parent_path/data/processed_promoters/$result/${FILENAME}_promoters.txt" &
done
)
