#!/bin/bash

# index for files
group=${1:-110}


PARENT_PATH=$(dirname $(greadlink -f $0))
result=${PARENT_PATH##*/}
echo $result
# Go back path
PARENT_PATH=${PARENT_PATH%/*}
PARENT_PATH=${PARENT_PATH%/*}
PARENT_PATH=${PARENT_PATH%/*}

# Find data directory
folder=$PARENT_PATH'/data/processed_sequencing/'$result

# Make directories for stored data
mkdir $PARENT_PATH'/data/processed_promoters/'
mkdir $PARENT_PATH'/data/processed_promoters/'$result

# Go to data
cd $folder

# Sequence file names
file_bc=$group'_R2.fastq.gz'
file_prom=$group'_R1.fastq.gz'

# Find barcodes
gunzip -c $file_bc | awk ' NR%4==2 {
        print $0;
    }
    NR%4==0 {
        print $0;
    }
    ' > $group'_barcodes.txt'

# Find promoters
gunzip -c $file_prom | awk ' NR%4==2 {
        print $0;
    }
    NR%4==0 {
        print $0;
    }
    ' > $group'_promoters.txt'

# Combine promoters and barcodes
paste $group'_barcodes.txt' $group'_promoters.txt' > $group'_combined.txt'

# Sort and count unique combinations
cat $group'_combined.txt' | awk 'NR%2==1 {print $0}' | sort | uniq -c | sort -bgr > $group'_collapsed.txt'

# Do some formatting
#tr -s "\t" " " < $group'_collapsed.txt' > $group'_collapsed.txt'

# Remove temporary files
rm $group'_promoters.txt'
rm $group'_barcodes.txt'
rm $group'_combined.txt'

awk '{
    print ">"$2"_"$1"\n"$3;
    print "";
}'  $group'_collapsed.txt' > $group'_collapsed.fasta'

mv $group'_collapsed.fasta' $parent_path'/data/processed_promoters/'$result/$group'_collapsed.fasta'