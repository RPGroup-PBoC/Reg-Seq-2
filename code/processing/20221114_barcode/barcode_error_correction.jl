using CSV, DataFrames, Distances, Printf

# Set path
dir = @__DIR__
home_dir = joinpath(split(dir, "/")[1:end-3])

df = CSV.read("/$home_dir/data/extracted_barcodes/20221114_barcode/$(file)_collapsed.txt", DataFrame, header=["count", "barcode"], ignorerepeated=true, delim=" ")

df = df[[length(x) == 20 for x in df.barcode], :]

df_singles = df[df.count .<= 2, :]
df_mult = df[df.count .> 2, :]

for (i, bc) in enumerate(df_singles.barcode)
    ind = findfirst(x -> hamming(bc, x) <=1, df_mult.barcode)
    if ~isnothing(ind)
        df_mult[ind, :count] += 1
    end
    if i%100==0
        println(@sprintf "%.4f%%  of barcodes processed.." i/nrow(df_singles) *100)
    end
end