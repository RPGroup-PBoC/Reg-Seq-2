---
status: Accepted
---

# 2022-05-14 Mapping Sequencing run

## Purpose

In this sequencing run we map the mutated promoters to barcodes.

## Platform

NextSeq (Caltech Sequencing facility)

## Template

## Primers used

PCR1 Forward and Sequencing Read 1 Primer:

TTGGTCCAGTGCCATGTTATCCCTGAATCTAGT

PCR1 Reverse and Sequencing Read 2 Primer:

ATACCTGTAGCTAAATCCCACCCGATGCTCGAC

PCR2 Forward Primer:

AATGATACGGCGACCACCGAGATCTACACTTGGTCCAGTGCCATGTTATCC

PCR2 Reverse Primers:

100 - CAAGCAGAAGACGGCATACGAGATCAGTGTATACCTGTAGCTAAATCCCACCC

110 - CAAGCAGAAGACGGCATACGAGATAGCCATATACCTGTAGCTAAATCCCACCC

201 - CAAGCAGAAGACGGCATACGAGATGTACTGATACCTGTAGCTAAATCCCACCC

204 - CAAGCAGAAGACGGCATACGAGATATCTCGATACCTGTAGCTAAATCCCACCC

Sequencing Index Primer:

GTCGAGCATCGGGTGGGATTTAGCTACAGGTAT

## Sequencing kit

NextSeq P2 flowcell, PE reads, 170 cycles on read1 and 30 cycles on read2.

## Materials

| **id** | **barcode-sequence** | **description** |
| :--: | :--: | :--: |
| 100 | ACACTG | Entire library |
| 110 | ATGGCT | Library Subpool |
| 201 | CAGTAC | Library Subpool |
| 204 | CGAGAT | Library Subpool |

## Processing
Store the sequencing files in

```
project
│   README.md  
│
└───data
│   └───sequencing
|       └───20220514_mapping
│           │   100_S1_R1_001.fastq.gz
│           │   100_S1_R2_001.fastq.gz
│           │   110_S2_R1_001.fastq.gz
│           │   110_S2_R1_001.fastq.gz
│           │   201_S3_R1_001.fastq.gz
│           │   201_S3_R1_001.fastq.gz
│           │   204_S4_R1_001.fastq.gz
│           │   204_S4_R1_001.fastq.gz

```

The files are processed using `fastp`. Execute `processing_seq.jl` to filter the sequencing files. This can be done by either starting the Julia REPL and using 

```julia
include("path/to/processing.jl")
```

or by running in the command line 

```
julia path/to/processing.jl
```


Then, we extract the barcodes and promoters from the sequencing data. In the following, for every script 

```
chmod +x extract_promoters.sh
```

to make the file executable. Then, simply run

```
./extract_promoters.sh <group number>
```

where `<group number>` is replaced by the group number of genes that is supposed to be processed. Here it is one of 100, 110, 201 and 204.
The script creates files containing each barcode and promoter pair, as well as their counts. The results will be stored in a `.fastq` file, which will be used to map the sequences to promoters.

To identify which promoter a sequence belongs to, we first need to prepare the data. Therefore, run the script 

```
run_bbmap.sh <group number>
```

Make sure that `bbmap` is in the correct location. To run this script, the current working directory has to be the folder of this experiment. Finally, to finish the alignment, run 

```
extract_gene_names.sh <group number>
```

in this folder. The script produces a file for each gene in the experiment, which contains the promoter variants as well as the counts.

The files can be combined again by running `filter_mapping.sh`, which also filters out all combinations of promoter and barcode that have less than 3 counts.

```
filter_mapping <group number>
```