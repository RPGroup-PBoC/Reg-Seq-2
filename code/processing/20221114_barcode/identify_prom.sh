gc=$1
awk 'NR==FNR{a[$1]=$0;next} ($3 in a){b=$3;$3="";print $0 a[b]}'  "../../barcodes/20220514_mapping/mapped_barcodes.csv" "$(gc)_combined.txt" > "$(gc)_identified.txt"
echo -e "cDNA_count gDNA_count barcode name mapping_count promoter" | cat - "$(gc)_identified.txt" > "$(gc)_identified_.txt"