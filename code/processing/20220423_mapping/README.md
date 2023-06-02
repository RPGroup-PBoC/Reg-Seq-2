---
status: Rejected
---

# 2022-04-23 Sequencing run

## Purpose
Mapping run for 119 promoters.

## Platform
NextSeq (Caltech Sequencing facility)

## Sequencing kit
NextSeq P2 flowcell, PE reads, 30 cycles on read1 and 170 cycles on read2.

## Materials

| **id** | **barcode-sequence** | **description** |
| :--: | :--: | :--: |
| 100 | CAGTGT | Entire library |
| 110 | AGCCAT | Library Subpool |
| 201 | GTACTG | Library Subpool |
| 204 | ATCTCG | Library Subpool |

## Notes and Observations
This NextSeq run failed. Read1 for the barcode had great quality and gave us many reads. However, read2 for the barcode had very poor quality, which we identified was an issue with the custom primers on the NextSeq flow cell. This problem will be avoided in future runs by the sequencing center.