# Experimental Design
In this folder you find all files needed to design the experiment. All files are Julia scripts which can be 
run in the project environment. For the environment setup, please follow the instructions on the frontpage
of the repository.

## Sequence Design

### Import and Processing of Databases

To identify the regulatory architecture of a gene, use the transcription start site to identify the region
which most likely contains binding sites for transcription factors. [Ecocyc](https://ecocyc.org/) is a database
which contains an annotated genome for *E. coli* K12 MG1655, which is the model organism in this project.
In this database we look for annotated promoters for the genes of interest. The script `get_gene_information_ecocyc.jl`
searches for promoters for a list of genes and returns a list of promoters and their TSS. Another database
is [RegulonDB](https://regulondb.ccg.unam.mx/). While these databases are mostly overlapping, there can be instances where
a promoter can be found in one database but not the other, so we search this database as well by running the script
`get_gene_information_regulonDB.jl`. To combine the results of the two searches, run the script `database_processing.jl`.
Finally, we design the sequences as explained in the methods section (link) using the script `design_sequences.jl`.
This script will ask for a prompt to choose the restriction enzymes that are used. One should use restriction enzymes 
which have a recognition site of length 6, otherwise the sequences won't have the correct length. The prompt also shows
which enzymes are used as default, which can be chosen by simply pressing `enter`.

All scripts can be run together by running the script `run_scripts.jl`.