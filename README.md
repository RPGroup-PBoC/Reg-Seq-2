# 1000_genes_ecoli

This is the Github Repository for the whole genome Reg-Seq project.
The code in this project is written in Julia, which can be downloaded [here](https://julialang.org/downloads/).


## Setting up Computational Environment

### Julia

To run the code in this project, you need to activate the custom environment, which can be done by starting Julia in this project folder. This can either be done by adding the Julia path as an Environment variable, or by starting the executable Julia file and navigating into the project folder (using the shell mode by typing `;`). Instructions on how Julia can be started from the command line can be found [here](https://julialang.org/downloads/platform/). Once the correct folder is selected, the working environment can be set by typing in the Julia REPL

```julia
julia> ] activate .
```

Once the environment is activated, all necessary packages can be installed with

```julia
julia> ] instantiate
```

To run a script use

```julia
julia> include("path/to/script.jl")
```

### Conda

During processing of sequencing data, we use the software package [`fastp`](https://github.com/OpenGene/fastp). This software can be installed using `conda`, therefore we need to set up an
appropriate `conda` environment. We provide a suitable `conda` environment that was used to run the processing for this project. The environment can be installed by running

```
conda env create -f environment.yml
```

Subsequently, the environment needs to be activated in order to be used. This can be done by

```
conda activate wgregseq
```

### BBmap

During processing of sequencing data, we use [bbmap](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/installation-guide/). It can be downloaded and used straight away. BBmap requires a working Java installation on the machine.

## Repository structure

### code

Contains all code files used to design experiments, process and analyze data and create figures.

### src

Contains the custom `Julia` software module.

### notebooks

Contains explanatory notebooks that walk through certain steps in the experimental design and data processing.

### test

Contains tests for software module.

### data

Contains supplementary data files.