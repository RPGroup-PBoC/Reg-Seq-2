#!/bin/bash
group=$1


result=${PWD##*/}
parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

parent_path=${parent_path%/*}
parent_path=${parent_path%/*}
parent_path=${parent_path%/*}
folder=$parent_path'/data/processed_sequencing/'$result
mkdir $parent_path'/data/processed_promoters/'
mkdir $parent_path'/data/processed_promoters/'$result

cd $folder
<<comment
file_bc=$group'_R2.fastq.gz'
file_prom=$group'_R1.fastq.gz'

gunzip -c $file_bc | awk ' NR%4==2 {
        print $0;
    }
    ' > $group'_barcodes.txt'

gunzip -c $file_prom | awk ' NR%4==2 {
        print $0;
    }
    ' > $group'_promoters.txt'
comment

paste $group'_barcodes.txt' $group'_promoters.txt' > $group'_combined.txt'
sort $group'_combined.txt' | uniq -c | sort -bgr > $group'_collapsed.txt'


rm $group'_promoters.txt'
rm $group'_barcodes.txt'
rm $group'_combined.txt'




