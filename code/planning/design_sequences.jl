using wgregseq, FASTX, DataFrames, CSV, BioSequences

# Set path
dir = @__DIR__
homedir = joinpath(split(dir, "/")[1:end-2])

# Import genome
r = open(FASTA.Reader, "/$homedir/data/ecocyc/mg1655_genome.fasta")
wt_sequence = [sequence(record) for record in r][1]

# Import gene list to generate sequences for
gene_table = CSV.read("/$homedir/data/100_genes.csv", DataFrame)

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

##
promoter_list[map(x-> "rspA" in x, promoter_list.genes), :]
##

# Get promoters for genes
df_list = DataFrame[]
for gene in gene_table.name
    _df = promoter_list[map(x -> gene in x, promoter_list.genes), [:promoter_ID, :tss, :direction, :gene_position, :genes]]
    push!(df_list, _df)
end
df = vcat(df_list...)

##
# Separate into promoters with or without tss
df_found_tss = df[.~isnan.(df.tss), :]
df_missing_tss = unique(df[isnan.(df.tss), :])
#df_missing_tss = df_missing_tss[map(x -> ~(x in df_found_tss.genes), df_missing_tss.genes), : ]


## find unique tss
gdf = groupby(df_found_tss, :tss)

# collapse genes with same tss
df_collapse = DataFrame()
for group in gdf
    _df = DataFrame(
        tss=unique(group.tss), 
        direction=unique(group.direction),
        genes=[unique(group.gene)],
        promoter_IDs=[unique(group.promoter_ID)],
    )
    df_collapse = vcat(df_collapse, _df)
end

df_collapse.tss = Int64.(df_collapse.tss)
df_collapse
CSV.write("/$homedir/data/gene_tss.csv", df_collapse)

## Design Sequences for genes with TSS
df_sequences = DataFrame()
for row in eachrow(df_collapse)
    tss = row.tss
    direction = row.direction
    genes = row.genes
    promoter_IDs = row.promoter_IDs
    seq = wgregseq.design.find_seq(tss, direction, 115, 45, wt_sequence)[1]
    mut_list = wgregseq.design.mutations_rand(seq, 0.1, 1500)
    _df = DataFrame(sequence=mut_list, genes=fill(genes, 1501), promoter_IDs=fill(promoter_IDs, 1501))
    df_sequences = vcat(df_sequences, _df)
end
df_sequences

## Create sequences for missing TSS

df_missing_tss
##
gdf = groupby(df_sequences, :promoter_IDs)
df_stack = DataFrame()

enzymes = ["SalI", "SbfI", "SacI", "NheI", "XbaI", "SpeI", "XhoI", "EcoRI", "ApaI", "ScaI", "NcoI", "MluI", "EcoRV", "BbsI", "BamHI", "AgeI"]

for group in gdf
    group[:, "sequence"] = wgregseq.design.add_primer(convert(Vector{BioSequences.LongDNASeq}, (group.sequence)), 100, "both")
    df_restriction = wgregseq.design.find_restriction_sites(enzymes, group[:, "sequence"])
    sort!(df_restriction, "sites")
    dict = Dict{Any, Any}(df_restriction.enzyme .=> df_restriction.sites)
    dict["gene"] = [unique(group.genes)[1]]
    dict["promoter"] = [unique(group.promoter_IDs)[1]]
    append!(df_stack, DataFrame(dict))
end


dict = Dict{Any, Any}(enzymes .=> sum.(eachcol(df_stack[!, enzymes])))
dict["gene"] = [String31["all"]]
dict["promoter"] = [String31["all"]]
append!(df_stack, DataFrame(dict))
println(df_stack)

##
