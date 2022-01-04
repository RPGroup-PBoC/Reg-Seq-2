using wgregseq, FASTX, DataFrames, CSV, BioSequences, CairoMakie, Jedi, StatsBase

Jedi.styles.default_makie!()

# Set path
dir = @__DIR__
homedir = joinpath(split(dir, "/")[1:end-2])

# Import promoter list and infer types that can be infered automatically
promoter_list = CSV.read(
    "/$homedir/data/promoter_list.csv", 
    DataFrame, 
    types=Dict(
        "TU_ID"=>String,
        "promoter_ID"=>String,
        "promoter"=>String,
        "tss"=>Float64,
        "direction"=>String
    )
)

# Define custom function for nice imports
Base.parse(::Type{Vector{String}}, x::String) = Vector{String}(filter(x-> x != ", ", split(x, "\""))[2:end-1])
function Base.parse(::Type{Vector{Float64}}, x::String)
    number = split(split(x, "[")[end][1:end-1], ", ")
    number_list = Float64[]
    for num in number
        if num != ""
            push!(number_list, parse(Float64, num))
        else
            return push!(number_list, NaN)
        end
    end
    return number_list

end
Base.parse(::Type{Vector{String}}, x::Missing) = String[]
Base.parse(::Type{Vector{Float64}}, x::Missing) = Float64[]

# Replace columns by nicer types
promoter_list.genes = parse.(Vector{String}, promoter_list.genes)
promoter_list.evidence = parse.(Vector{String}, promoter_list.evidence)
promoter_list.gene_position = parse.(Vector{Float64}, promoter_list.gene_position)

promoter_list

##
found_tss = promoter_list[(promoter_list.TU_ID .!= "None") .& .~(isnan.(promoter_list.tss)), :]
no_tss = promoter_list[(promoter_list.TU_ID .!= "None") .& isnan.(promoter_list.tss), :]


found_tss

min_distance = Float64[]
for i in 1:nrow(found_tss)
    min_dist = minimum(abs.(found_tss.gene_position[i] .- found_tss.tss[i]))
    push!(min_distance, min_dist)
end 

found_tss[!, "distance"] = min_distance
sort!(found_tss, "distance", rev=true)
println(first(found_tss[.~isnan.(found_tss.distance), ["genes", "distance"]], 15))
##
fig = Figure(resolution=(300, 300))
ax = Axis(fig[1, 1], xscale=log10)
ax.xlabel = "Bases"
ax.ylabel = "ECDF"
ax.title = "Distance of TSS to coding region"
ind = findfirst(x -> x > 0, sort(min_distance))

lines!(
    ax, 
    sort(min_distance)[ind:end], 
    collect(ind:length(min_distance)) ./ length(min_distance),
    linewidth=3
    
)

fig
save("figures/tss_distance.pdf", fig)

