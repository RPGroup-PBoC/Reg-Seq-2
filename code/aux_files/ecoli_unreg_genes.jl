using DataFrames, wgregseq, CSV, CairoMakie

# Set path
dir = @__DIR__
home_dir = joinpath(split(dir, "/")[1:end-2])

# Set plotting style
wgregseq.plotting_style.default_makie!()

# Import genes
promoter_list = CSV.read(
    "/$home_dir/data/promoter_list_processed.csv", 
    DataFrames.DataFrame, 
    types=Dict(
        "promoter"=>String,
        "tss"=>Float64,
        "direction"=>String
    )
)

# Replace columns by nicer types
promoter_list.genes = parse.(Vector{String}, promoter_list.genes)
promoter_list.gene_position = parse.(Vector{Float64}, promoter_list.gene_position)
promoter_list.evidence = parse.(Vector{String}, promoter_list.evidence)
promoter_list.TF_IDS = parse.(Vector{String}, promoter_list.TF_IDS)

# Find all unique gene names
all_genes = all_gene_list.gene |> unique

# Find all genes that have a promoter annotated 
genes_w_promoter = vcat(promoter_list.genes...) |> unique 
println("Genes with promoter in ecocyc $(length(genes_w_promoter) / length(all_genes))")


# Find all promoters with transcriptional annotation and get genes regulated by promoter
promoter_w_tfs = promoter_list[map(x -> x[1] != "none", promoter_list.TF_IDS), :genes]
# Find unique genes names 
genes_w_tf = vcat(promoter_w_tfs...) |> unique
println("Genes with TF annotated in ecocyc $(length(genes_w_tf) / length(all_genes))")