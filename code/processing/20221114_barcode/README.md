---
status: Accepted
---

# 2022-11-14 Barcode Sequencing run

## Purpose

In this sequencing run we count barcodes for a library of 119 promoters that have been grown in four growth conditions, mentioned below.

## Platform

NovaSeq (Stanford Sequencing facility)

## Template

## Primers used

Sequencing Read 1 Primer:

ATACCTGTAGCTAAATCCCACCCGATGCTCGAC


Sequencing Index Primer:

GTCGAGCATCGGGTGGGATTTAGCTACAGGTAT

## Sequencing kit

NovaSeq S2 flowcell, SE reads, 26 cycles on read1

## Materials

| **id** | **primer index** | **barcode-sequence** | **description** |
| :--: | :--: | :--: | :--:|
LB_cDNA	|	SC94 |	GCATAG | RNA-Seq in LB
Etoh_cDNA	|	SC93 |	CATTCG | RNA-Seq in LB + 5% Ethanol
Gluc_cDNA	|	SC95 |	ACTAGC | RNA-Seq in M9 minimal media + 0.5% Glucose
Xy_Ar_cDNA	|	SC96 |	CAGTAC | RNA-Seq in M9 minimal media + 0.5% Xylose + 0.5% Arabinose
LB_gDNA	|	SC90 |	ATGGCT | gDNA-Seq in LB
Etoh_gDNA	|	SC89 |	CACAGT | gDNA-Seq in LB + 5% Ethanol
Gluc_gDNA	|	SC91 |	CGAGAT | gDNA-Seq in M9 minimal media + 0.5% Glucose
Xy_Ar_gDNA	|	SC92 |	ACACTG | gDNA-Seq in M9 minimal media + 0.5% Xylose + 0.5% Arabinose

## Processing
Store the sequencing files in the following format:

```
project
│   README.md  
│
└───data
│   └───sequencing
|       └───20221114_barcode
│           │   LB_cDNA_S1_R1_001.fastq.gz
│           │   XyAr_gDNA_S8_R1_001.fastq.gz
│           │   LB_gDNA_S5_R1_001.fastq.gz
│           │   Etoh_cDNA_S2_R1_001.fastq.gz
│           │   XyAr_cDNA_S4_R1_001.fastq.gz
│           │   Etoh_gDNA_S6_R1_001.fastq.gz
│           │   Gluc_cDNA_S3_R1_001.fastq.gz
│           │   Gluc_gDNA_S7_R1_001.fastq.gz

```

The files are processed using `fastp`. The first step is to make the bash scripts in this folder executable. This can be done by using 

```
chmod +x *.sh
```

Next, we filter the sequences for quality scores and trim the trailing six bases from each read, since the first 20 bases of the read are the barcode, and the following bases are conserved. This is done by running the `processing_seq.sh` script in this folder. Run 

```
./processing_seq.sh
```

Then, we extract the barcodes and promoters from the sequencing data. Run

```
./extract_barcodes.sh
```

The script creates files containing each barcode, as well as their counts. The results will be stored in a list (and `.fasta` file), which will be used to map the sequences to promoters.
**Insert here part about barcode qc**