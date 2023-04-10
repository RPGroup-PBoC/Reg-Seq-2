# Processing of Mapping Dataset
The raw data is stored as `fastq` files for read 1 and read 2. Read 1 contains the promoter variant, while read 2 contains the associated random barcode. Instructions on how to install the software necessary to process the data can be found on the [main page](https://github.com/RPGroup-PBoC/1000_genes_ecoli) of this repository.

## Quality Filtering
First, the data is filtered on the read quality. Requires `fastp` to be installed. `1_quality_filter.sh`

## Extraction of Promoter and Barcode sequence
For each read, promoter sequence and barcode sequence are extracted from the `fastq` files and stored in `fasta` file that can be used to identify the promoter each variant is derived from. Requires `GNU parallel` to be installed. `2_extract_promoters.sh`

## Identification of promoter variants
Identify promoters. Requires `bbmap` to be installed. `3_extract_gene_names.sh`.