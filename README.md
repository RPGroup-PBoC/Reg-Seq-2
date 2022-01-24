# 1000_genes_ecoli

This is the Github Repository for the whole genome Reg-Seq project.
The code in this project is written in Julia, which can be downloaded [here](https://julialang.org/downloads/).

To run the code in this project, you need to activate the custom environment, which can be done by starting Julia in this project folder and typing in the Julia REPL

```julia
julia> ] activate .
```

To run a script use 

```julia
julia> include("path/to/script.jl")
```

## code

Contains scripts that can be run to reproduce the experiment, from designing sequences to analyzing sequencing results and making figures.

## data

Contains data files that are either needed for preparing the experiment or is the destination for data files from experiments (not included
in this repository, can be found at (link)).

## src

Contains software module in Julia used in this project.

## test

Should contain tests of the software module.