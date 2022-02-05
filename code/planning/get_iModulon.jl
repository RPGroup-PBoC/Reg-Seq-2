using CSV, DataFrames, wgregseq, CairoMakie

wgregseq.plotting_style.default_makie!()
##
# Open Gene information
dir = @__DIR__
home_dir = joinpath(split(dir, '/')[1:end-2]...)


## Metadata for samples
df_samples = CSV.read("/$home_dir/data/iModulon/sample_table.csv", DataFrames.DataFrame)
names(df_samples)

# Fill missing values
df_samples[!, ["Carbon Source (g/L)"]] = coalesce.(df_samples[!, ["Carbon Source (g/L)"]], "none")
df_samples[!, ["Strain Description"]] = coalesce.(df_samples[!, ["Strain Description"]], "none")

# Take only MG1655
df_samples_MG = df_samples[
    Array(df_samples[!, "Strain Description"]) .== "Escherichia coli K-12 MG1655",
     :]

sample_IDs = df_samples_MG.Column1
df_samples[!, "Strain Description"] |> unique
df_samples_MG[!, ["Carbon Source (g/L)"]] |> unique |> println
df_samples_MG[!, ["Nitrogen Source (g/L)"]] |> unique  |> println
df_samples_MG[!, ["Electron Acceptor"]] |> unique |> println
df_samples_MG[!, ["Trace Element Mixture"]] |> unique |> println
## iModulon metadata
df_im = CSV.read("/$home_dir/data/iModulon/iM_table.csv", DataFrames.DataFrame)
names(df_im)
df_im.category |> unique
df_im.n_genes |> sum
df_im[43, :]
println(df_im[df_im.category .== "Carbon Metabolism", ["name", "regulator_readable", "n_genes", ]])

println(filter(x-> occursin("lac", x), df_im.name))

## Gene - iModulon list
df_gene_im_list = CSV.read("/$home_dir/data/iModulon/gene_presence_list.csv", DataFrame)
df_gene_im_list.Gene |> unique |> length
## Gene - iModulon boolean table
df_gene_im_table = CSV.read("/$home_dir/data/iModulon/gene_presence_matrix.csv", DataFrame)

## Raw Counts
df_counts = CSV.read("/$home_dir/data/iModulon/Precise2dot0counts.csv", DataFrame)

## Normalized reads 
df_tpm = CSV.read("/$home_dir/data/iModulon/log_tpm.csv", DataFrame)
df_tpm = df_tpm[!, vcat(["Column1"], sample_IDs)]

##
df_M = CSV.read("/$home_dir/data/iModulon/M.csv", DataFrame)

# Import gene information
df_genes = CSV.read("/$home_dir/data/all_genes_table.csv", DataFrame)


## Test lacI 
lacZ_acc = df_genes[df_genes.gene .== "lacZ", "accession"][1]
lacZ_iM = df_gene_im_list[df_gene_im_list.Gene .== lacZ_acc, "iModulon"][1]
lacZ_iM_genes = df_gene_im_list[df_gene_im_list.iModulon .== lacZ_iM, "Gene"]
df_genes[map(x -> x in lacZ_iM_genes, df_genes.accession), "gene"]

## Normalized reads 
df = CSV.read("/$home_dir/data/iModulon/A.csv", DataFrames.DataFrame)
Array(df_samples[!, "Base Media"]) .== "M9"
df_samples = df_samples[
    (Array(df_samples[!, "Strain Description"]) .== "Escherichia coli K-12 MG1655") .&
    (Array(df_samples[!, "Base Media"]) .== "M9"), :]
df = df[!, df_samples.Column1]

x = df[24, names(df)[2:end]]
gal_samples = df_samples[df_samples[!, "Carbon Source (g/L)"] .== "galactose(4)", "Column1"]

sig_ind = Array(x) .> 2.5
names(x[sig_ind]) 
df_samples.Column1
_df = DataFrame()
for y in names(x[sig_ind])

    append!(_df, df_samples[(df_samples.Column1 .== y), :])
end
_df



fig = Figure()
ax = Axis(fig[1, 1])
scatter!(Array(x))
fig

## Knockouts

gdf = groupby(df_im, "category")
gene_cat_df = DataFrame()
for _df in gdf
    gene_list = []
    for i in _df.k
        genes = df_gene_im_list[df_gene_im_list.iModulon .== i, "Gene"]
        gene_names = df_genes[map(x -> x in genes, df_genes.accession), "gene"]
        push!(gene_list, gene_names...)
    end
    append!(gene_cat_df, DataFrame(genes=gene_list, category=fill(unique(_df.category)[1], length(gene_list))))
end
gene_cat_df[gene_cat_df.genes .== "bhsA", :]

gene_cat_df

cat_per_gene = sort(combine(groupby(unique(gene_cat_df), "genes"), nrow), "nrow", rev=true)
mean(cat_per_gene.nrow)

gene_cat_df[map(x -> occursin("lac", x), gene_cat_df.genes), :]


##
unique_genes = df_gene_im_list.Gene |> unique 
genes_precise = df_genes[map(x -> x in unique_genes, df_genes.accession), "gene"]
df_heinemann = DataFrames.DataFrame(Pandas.read_excel("data/heinemann_data.xlsx", sheet_name="Table S6", header=[2]))
genes_heinemann = df_heinemann.Gene |> unique

sum(map(x -> x in genes_heinemann, genes_precise))
# compare synomyms

## Some genes may have the wrong synomym
all_gene_list = CSV.read(
    "/$home_dir/data/all_genes_table.csv", 
    DataFrame, 
    types=Dict(
        "ID"=>String,
        "gene"=>String,
        "gene_position"=>Float64,
        "direction"=>String
    )
)

all_gene_list.synonyms = parse.(Vector{String}, all_gene_list.synonyms)
all_gene_list.transcription_units = parse.(Vector{String}, all_gene_list.transcription_units)

for l in [genes_precise, genes_heinemann]
    for i in 1:length(l)
        name = l[i]
        if name ∉ all_gene_list.gene
            if length(all_gene_list[map(x -> name in x, all_gene_list.synonyms), :gene]) == 0
                println(name)
            else
                syn = all_gene_list[map(x -> name in x, all_gene_list.synonyms), :gene][1]
                l[i] = syn
            end
        end
    end
end

sum(map(x -> x in genes_heinemann, genes_precise))