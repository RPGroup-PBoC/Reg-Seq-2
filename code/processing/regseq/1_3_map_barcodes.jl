using DataFrames, CSV, BioSequences, wgregseq, CairoMakie, Distances, Statistics


# Set path
dir = @__DIR__
home_dir = joinpath(split(dir, "/")[1:end-3])

df_wt = CSV.read("/$home_dir/data/wtsequences_regseq.csv", DataFrame)

df_out = CSV.read(
    "/Volumes/rp_lab_ext/og_regseq_data/mapping_counted.csv", 
    DataFrame, 
    ignorerepeated=true,
    delim="\t",
    header=["counts", "promoter", "barcode"]
    )



names = String[]
nmuts = Int[]
promoters = String[]
barcodes = String[]
counts = Int[]
gdf = groupby(df_out, :barcode)
for _df in gdf
    if nrow(_df) > 1
        skip
    end
    x = LongDNA{4}(_df.promoter[1])
    distances = [hamming(x, y) for y in LongDNA{4}.(df_wt.geneseq)]
    closest_seq =  argmin(distances)
    push!(names, df_wt.name[closest_seq])
    push!(nmuts, minimum(distances))
    push!(promoters, _df.promoter[1])
    push!(barcodes, _df.barcode[1])
    push!(counts, _df.counts[1])
end


df_out = DataFrame(barcode=barcodes, promoter=promoters, count=counts, name=names, nmut=nmuts)
CSV.write("/Volumes/rp_lab_ext/og_regseq_data/mapping_identified.csv", df_out)
display(first(df_out, 50))

for gene in df_out.name |> unique
    CSV.write(
        "/Volumes/rp_lab_ext/og_regseq_data/promoter_bc_key_per_gene/$(gene)_identified.csv", 
        df_out[
            (df_out.name .== gene) .& (df_out.count .> 1), :])
end