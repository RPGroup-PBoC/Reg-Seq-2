using wgregseq, FASTX, DataFrames, CSV, BioSequences, CairoMakie,  StatsBase

# Set path
dir = @__DIR__
home_dir = joinpath(split(dir, "/")[1:end-3])

df = CSV.read("/$home_dir/data/twist_orders/TWIST_sequences_niko_30000_new.csv", DataFrame)
insertcols!(df, 1, :promoter_name => map(x -> split(x, '_')[1], df.name))

df_negs = df[df.promoter_name .== "neg", :]
insertcols!(df_negs, 1, :promoter_seq => map(x -> x[21:180], df_negs.sequence))

df = df[df.promoter_name .!= "neg", :]

println(df.promoter_name |> unique)