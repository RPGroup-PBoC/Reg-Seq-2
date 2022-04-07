N=${1:-1}
FILES=$(find . -name '*.fastq.gz')
(
for file in $FILES
do
  ((i=i%N)); ((i++==0)) && wait
  FILENAME=${file%.fastq.gz}
  FILENAME=${FILENAME#*/}
  gunzip -c $file | awk /CCGCATTG/ | awk 'NR%4==0 { if ((length($0) == 142) || length($0) == 150) { bc = substr($0,0,20); prom=substr($0,length($0)-10,10); print bc "," prom } }' > "${FILENAME}_barcodes.txt" &
done
)
