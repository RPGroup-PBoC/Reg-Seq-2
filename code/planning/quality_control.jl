using wgregseq, DataFrames, CSV, BioSequences, CairoMakie, Statistics, StatsBase

wgregseq.plotting_style.default_makie!()

# Set path
dir = @__DIR__
home_dir = joinpath(split(dir, "/")[1:end-2])

##
filename = "2022-02-15_twist_order.csv"

df = CSV.read(
    "/$home_dir/data/twist_orders/$filename",
    types=Dict(
        "fwd_primer" => String,
        "rev_primer1" => String,
        "rev_primer2" => String,
        "rev_primer3" => String,
    ),
    DataFrame
)
df.fwd_primer = parse.(Tuple{Int64, Tuple{Int64, Int64}}, df.fwd_primer)
df.rev_primer1 = parse.(Tuple{Int64, Tuple{Int64, Int64}}, df.rev_primer1)
df.rev_primer2 = parse.(Tuple{Int64, Tuple{Int64, Int64}}, df.rev_primer2)
df.rev_primer3 = parse.(Tuple{Int64, Tuple{Int64, Int64}}, df.rev_primer3)


df.genes = parse.(Vector{String}, df.genes)
df.sequence = [LongDNA{4}(seq) for seq in df.sequence]
##
#wgregseq.quality_control.check_dataframe(df, site_start=27, site_end=186, print_results=false)
wgregseq.quality_control.check_primers_off_target(df)

##
gdf = groupby(df, :rev_primer3)
for _df in gdf
    println([pr[1] for pr in _df.rev_primer3 |> unique])
    println([pr[1] for pr in _df.rev_primer2 |> unique])
    println(nrow(_df))
    println(_df.promoter |> unique |> length)
    println()
end