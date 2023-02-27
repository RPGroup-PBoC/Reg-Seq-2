using CSV, DataFrames, Statistics

##
# Open Gene information
dir = @__DIR__
home_dir = joinpath(split(dir, '/')[1:end-3]...)

file = split(dir, '/')[end]
begin
    df = CSV.read("/$home_dir/data/barcodes/$file/mapping_identified.csv", DataFrame, header=["promoter", "barcode", "count", "name"])
    dropmissing!(df)

    # counts per promoter
    println("counts per promoter")
    df_counts = combine(groupby(df, :name), :count => sum)
    sort!(df_counts, :count_sum, rev=true)
    display(df_counts)
    println()

    df = df[df.count .> 1, :]
    cdf = combine(groupby(df, :barcode), nrow)
    bc_count = Dict(cdf.barcode .=> cdf.nrow)
    df = df[[bc_count[x] == 1 for x in df.barcode], :]
    

    num_unique(x) = unique(x) |> length
    # number of unique promoters found
    println("number of unique promoters found")
    df_counts = combine(groupby(df, :name), :promoter => num_unique)
    sort!(df_counts, :promoter_num_unique, rev=true)
    display(df_counts)
    println()

    # average number of barcodes per promoter
    println("average number of barcodes per promoter")
    df_counts = combine(groupby(df, [:name, :promoter]), nrow)
    df_counts = combine(groupby(df_counts, :name), :nrow => mean => :mean_barcodes)
    sort!(df_counts, :mean_barcodes, rev=true)
    display(df_counts)
    println()

    # unique promoters with more than or equal to 10 barcodes
    println("unique promoters with more than or equal to 10 barcodes")
    df_counts = combine(groupby(df, [:name, :promoter]), nrow)
    df_counts = df_counts[df_counts.nrow .> 9, :]
    df_counts = combine(groupby(df_counts, :name), :promoter => num_unique)
    sort!(df_counts, :promoter_num_unique, rev=true)
    display(df_counts)

end





