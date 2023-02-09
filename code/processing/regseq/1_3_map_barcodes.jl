using DataFrames, CSV, BioSequences, wgregseq, CairoMakie, Distances, Statistics


# Set path
dir = @__DIR__
home_dir = joinpath(split(dir, "/")[1:end-3])

df_wt = CSV.read("/$home_dir/data/wtsequences_regseq.csv", DataFrame)
df_groups = CSV.read("/$home_dir/data/genetogroupnum", DataFrame, select=[2, 3])
gene_group_dict = Dict(df_groups.genename .=> df_groups.pnum)

for file in readdir("/home/tom/git/1000_genes_ecoli/data/mapping_sequences"; join=true)
    println(file)
    if ~occursin("mapping_counted.csv", file)
        continue
    end

    df_in = CSV.read(
        file, 
        DataFrame, 
        ignorerepeated=true,
        delim="\t",
        header=["counts", "promoter", "barcode"]
        )
    group = split(split(file, '/')[end], '_')[1]
    names = String[]
    nmuts = Int[]
    promoters = String[]
    barcodes = String[]
    counts = Int[]
    gdf = groupby(df_in, :barcode)
    for _df in gdf
        x = LongDNA{4}(_df.promoter[1])
        distances = [hamming(x, y) for y in LongDNA{4}.(df_wt.geneseq)]
        closest_seq =  argmin(distances)
        if string(gene_group_dict[split(df_wt.name[closest_seq], '_')[1]]) == group
            push!(names, df_wt.name[closest_seq])
            push!(nmuts, minimum(distances))
            push!(promoters, _df.promoter[1])
            push!(barcodes, _df.barcode[1])
            push!(counts, _df.counts[1])
        end
    end


    df_out = DataFrame(barcode=barcodes, promoter=promoters, count=counts, name=names, nmut=nmuts)
    CSV.write("/home/tom/git/1000_genes_ecoli/data/mapping_sequences/$(group)_mapping_identified.csv", df_out, delim=" ")
    display(first(df_out, 50))

    for gene in df_out.name |> unique
        CSV.write(
            "/home/tom/git/1000_genes_ecoli/data/promoter_bc_key_per_gene/$(gene)_identified.csv", 
            df_out[
                (df_out.name .== gene) .& (df_out.count .> 1), :])
    end
end