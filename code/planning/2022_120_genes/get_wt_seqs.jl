using CSV, DataFrames, wgregseq

# Set path
dir = @__DIR__
home_dir = joinpath(split(dir, "/")[1:end-3])

df = wgregseq.utils.import_twist_order(
    "/$home_dir/data/twist_orders/2022-02-15_twist_order.csv"
)

gdf = groupby(df, :promoter)

for _df in gdf
    println(">", _df[1, :promoter])
    println(_df[1, :sequence][27:186])
    println()
end

names(df)