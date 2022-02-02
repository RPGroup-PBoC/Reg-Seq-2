using wgregseq, FASTX, DataFrames, CSV, BioSequences, CairoMakie,  StatsBase

# Set path
dir = @__DIR__
home_dir = joinpath(split(dir, "/")[1:end-2])

# Import promoter list and infer types that can be infered automatically

# Promoters from Ecocyc
promoter_list_ecocyc = CSV.read(
    "/$home_dir/data/promoter_list_ecocyc.csv", 
    DataFrames.DataFrame, 
    types=Dict(
        "promoter"=>String,
        "tss"=>Float64,
        "direction"=>String
    )
)

# Promoters from RegulonDB
promoter_list_regulonDB = CSV.read(
    "/$home_dir/data/promoter_list_regulon_DB.csv", 
    DataFrames.DataFrame, 
    types=Dict(
        "promoter"=>String,
        "tss"=>Float64,
        "direction"=>String
    )
)


# Replace columns by nicer types
promoter_list_ecocyc.genes = parse.(Vector{String}, promoter_list_ecocyc.genes)
promoter_list_ecocyc.gene_position = parse.(Vector{Float64}, promoter_list_ecocyc.gene_position)

promoter_list_regulonDB.genes = parse.(Vector{String}, promoter_list_regulonDB.genes)
promoter_list_regulonDB.gene_position = parse.(Vector{Float64}, promoter_list_regulonDB.gene_position)


# Join the datasets
df = promoter_list_ecocyc[(.~ isnan.(promoter_list_ecocyc.tss)) .& (map(x-> x != ["None"], promoter_list_ecocyc.genes)), :]
promoter_list_regulonDB.tss = coalesce.(promoter_list_regulonDB.tss, NaN)
append!(df, promoter_list_regulonDB)

##

# Drop missing entries
dropmissing!(df)
df.direction = convert(Vector{String}, df.direction)

# Sort gene entries to remove duplicates
for i in 1:nrow(df)
    genes, positions = df[i, ["genes", "gene_position"]]
    index = sortperm(genes)
    df[i, ["genes", "gene_position"]] = genes[index], positions[index]
end

# Alphabetically sort dataframe and remove duplicates
sort!(df, "genes")
df = unique(df)

# Last 5 entries are artifacts
df = df[5:end, :]

# Store list
CSV.write("/$home_dir/data/promoter_list_processed.csv", df)

##
df_genes = CSV.read("/$home_dir/data/all_genes_table.csv", DataFrames.DataFrame)[1:end-1, :]
# Import promoters with genes annotated
df = df[map(x -> length(x)>0, df.genes), :]

# Find genes which have a promoter annotated and which don't
gdf = groupby(df, "promoter")
gene_promoter = []
for group in gdf
    push!(gene_promoter, group.genes...)
end
genes_w_promoter = unique(vcat(gene_promoter...))
genes_wo_promoter = filter(x -> x ∉ genes_w_promoter, df_genes.gene)
println("Percentage of genes without a promoter: $(length(genes_wo_promoter)/(length(genes_wo_promoter)+length(genes_w_promoter)))")