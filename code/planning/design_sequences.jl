using wgregseq, FASTX, DataFrames, CSV, BioSequences, CairoMakie, Statistics, StatsBase

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
    "/$homedir/data/promoter_list_processed.csv", 
    DataFrame, 
    types=Dict(
        "promoter"=>String,
        "tss"=>Float64,
        "direction"=>String
    )
)

operons_without_promoters = CSV.read(
    "/$homedir/data/operons_without_promoters.csv", 
    DataFrame, 
    types=Dict(
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
promoter_list.gene_position = parse.(Vector{Float64}, promoter_list.gene_position)

operons_without_promoters.genes = parse.(Vector{String}, operons_without_promoters.genes)
operons_without_promoters.gene_position = parse.(Vector{Float64}, operons_without_promoters.gene_position)

promoter_list
operons_without_promoters
## Some genes may have the wrong synomym
all_gene_list = CSV.read(
    "/$homedir/data/all_genes_table.csv", 
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

for i in 1:nrow(gene_table)
    name = String(gene_table[i, "name"])
    if name ∉ all_gene_list.gene
        syn = all_gene_list[map(x -> name in x, all_gene_list.synonyms), :gene][1]
        gene_table[i, "name"] = syn
    end
end


all_gene_list

println("Data imported.")
## Get promoters for genes
df_list = DataFrame[]
genes_no_tss = []
for gene in gene_table.name
    _df = promoter_list[map(x -> gene in x, promoter_list.genes), [:tss, :direction, :gene_position, :genes, :promoter]]
    if nrow(_df) > 0
        push!(df_list, _df)
    else
        push!(genes_no_tss, gene)
    end
end
df = vcat(df_list...) |> unique
println("Promoters found for genes.")
println()

## For genes with multiple promoters, take the strongest predicted one
function find_best_promoter(df)
    min_tss = df.tss |> minimum
    max_tss = df.tss |> maximum
    sequence = wt_sequence[Int(min_tss)-115:Int(max_tss)+115]
    
    p = wgregseq.promoter_finder.Promoter_Calculator()
    r = p(sequence)
    if df.direction[1] == "-"
        sites = [r["Reverse_Predictions_per_TSS"][x - (Int(min_tss)-115)] for x in df.tss]
    else
        sites = [r["Forward_Predictions_per_TSS"][x - (Int(min_tss)-115)] for x in df.tss]
    end
    ind = argmin([site["dG_total"] for site in sites])
    return df[ind, :]
end

_df = DataFrame()
for gene_list in df.genes |> unique
    temp_df = df[map(x -> x == gene_list, df.genes), :]
    if nrow(temp_df) > 1
        append!(_df, DataFrame(find_best_promoter(temp_df)))
    else
        append!(_df, temp_df)
    end
end

df = _df

println("Best promoter found for genes with multiple promoters.")
println()
## Missing TSS
# Group genes into operons / transcription units
df_no_prom = DataFrame()
for gene in genes_no_tss
    append!(df_no_prom, operons_without_promoters[map(x -> gene in x, operons_without_promoters.genes), :])
end
unique!(df_no_prom)

# Look for these units in Urtecho et al. 2020 dataset

# import Urtecho data
urtecho_tss = CSV.read(
    "/$homedir/data/urtecho_2020/tss_operon_regulation.txt", 
    DataFrame 
)


function occursin_operon(gene, operon)
    split_operon = split(operon, "-")
    a, b = gene[1:3], gene[4]
    return prod(occursin.(a, split_operon) .* occursin.(b, split_operon))
end

# Find operon in Urtecho data
df_tss_urtecho = DataFrame()
# Make array to store indeces of genes to be removed from list
delete_index_list = Int64[]
for i in 1:nrow(df_no_prom)
    genes = df_no_prom[i, "genes"]
    gene_position = df_no_prom[i, "gene_position"]
    for gene in genes
        operons = filter(x -> prod(occursin_operon.(gene, x)), urtecho_tss.operon)
        if (operons |> unique |> length) > 1
            throw(ErrorException("More than one operon for genes: $(genes)"))
        elseif  (operons |> unique |> length) == 0
            println("No operons for genes $(genes) in Urtecho dataset.")
        else
            operon = unique(operons)[1]
            temp = urtecho_tss[(urtecho_tss.operon .== operon) .& ((urtecho_tss.active .== "active")), :]
            if nrow(temp) != 0
                insertcols!(temp, 2, :genes => fill(genes, nrow(temp)))
                insertcols!(temp, 2, :gene_position => fill(gene_position, nrow(temp)))
                rename!(temp, "tss_strand" => "direction", "tss_position"=> "tss", "tss_name"=>"promoter")
                append!(df_tss_urtecho, temp)
                push!(delete_index_list, i)
            end
        end
    end
end
append!(df, df_tss_urtecho[:, ["genes", "tss", "direction", "gene_position", "promoter"]])
df_no_prom = df_no_prom[Not(delete_index_list), :]
println("Took TSS from Urtecho data set.")
println()
## Now we tackle the genes for which we have not found a TSS yet

df_no_prom

p = wgregseq.promoter_finder.Promoter_Calculator()
tss_list = Int64[]
for i in 1:nrow(df_no_prom)
    if df_no_prom[i, "direction"] == "+"
        ind = Int64(minimum(df_no_prom[i, "gene_position"]))
        sequence = wt_sequence[ind-500:ind]
        r = p(sequence)["Forward_Predictions_per_TSS"]
        _x = [(key, r[key]["dG_total"]) for key in keys(r) |> collect]
        tss = _x[argmin([x[2] for x in _x])][1] + ind - 500
    else
        ind = Int64(maximum(df_no_prom[i, "gene_position"]))
        sequence = wt_sequence[ind:ind+500]
        r = p(sequence)["Reverse_Predictions_per_TSS"] 
        _x = [(key, r[key]["dG_total"]) for key in keys(r) |> collect]
        tss = _x[argmin([x[2] for x in _x])][1] + ind
    end
    push!(tss_list, tss)
end
tss_list
insertcols!(df_no_prom, 2, :tss =>tss_list)
insertcols!(df_no_prom, 5, :promoter =>fill("predicted", nrow(df_no_prom)))
append!(df, df_no_prom)

df

println("Predicted TSS using La Fleur model.")
println()
## Design Sequences for genes with TSS
df_sequences = DataFrame()
for row in eachrow(df)
    tss = Int64(row.tss)
    direction = row.direction
    genes = row.genes
    promoter = row.promoter
    seq = wgregseq.design.find_seq(tss, direction, 115, 45, wt_sequence)[1]
    mut_list = wgregseq.design.mutations_rand(seq, 0.1, 1500)
    _df = DataFrame(sequence=mut_list, genes=fill(genes, 1501), promoter=fill(promoter, 1501))
    df_sequences = vcat(df_sequences, _df)
end
df_sequences

if any(length.(df_sequences.sequence) .!= 160)
    throw(ErrorException("Not all sequences are 160bp!"))
else
    println("Mutated sequences created!")
    println()
end

## Check sequences for restriction sites 
gdf = groupby(deepcopy(df_sequences), :genes)
df_stack = DataFrame()

enzymes = ["SalI", "SbfI", "SacI", "NheI", "XbaI", "SpeI", "XhoI", "EcoRI", "ApaI", "ScaI", "NcoI", "MluI", "EcoRV", "BbsI", "BamHI", "AgeI"]

for group in gdf
    group[:, "sequence"] = wgregseq.design.add_primer(convert(Vector{BioSequences.LongDNASeq}, (group.sequence)), 100, "both")
    df_restriction = wgregseq.design.find_restriction_sites(enzymes, group[:, "sequence"])
    sort!(df_restriction, "sites")
    dict = Dict{Any, Any}(df_restriction.enzyme .=> df_restriction.sites)
    dict["gene"] = [unique(group.genes)[1]]
    dict["promoter"] = [unique(group.promoter)[1]]
    append!(df_stack, DataFrame(dict))
end
df_stack

dict = Dict{Any, Any}(enzymes .=> sum.(eachcol(df_stack[!, enzymes])))
dict["gene"] = [String31["all"]]
dict["promoter"] = "all"
append!(df_stack, DataFrame(dict))

## Add restriction sites
df_sequences.sequence = wgregseq.design.add_re_sites.(df_sequences.sequence, "SacI", "SalI")
if any(length.(df_sequences.sequence) .!= 172)
    throw(ErrorException("Not all sequences are 172 after adding primers bp!"))
end
## Add primers

df_sequences.sequence = wgregseq.design.add_primer(df_sequences.sequence, 100)

# Confirm that the primers are correct
fwd_primer = wgregseq.design.import_primer(100, "fwd")
rev_primer = wgregseq.design.import_primer(100, "rev")
if any([seq[1:20] != fwd_primer for seq in df_sequences.sequence])
    throw(ErrorException("Not all sequences have the right forward primer!"))
elseif any([seq[end-19:end] != rev_primer for seq in df_sequences.sequence])
    println("Forward primer and first reverse primer added.")
    println()
end

## Add reverse primers
n_per_group = 5
cop_df = deepcopy(df_sequences)
insertcols!(cop_df, 4, :rev_primer2 => fill(0, nrow(cop_df)))
insertcols!(cop_df, 5, :rev_primer3 => fill(0, nrow(cop_df)))

gdf = groupby(cop_df, :genes)

i = 1
while i * n_per_group < length(gdf)
    for _df in gdf[1+(i-1)*n_per_group:i*n_per_group]
        _df[:, "sequence"] = wgregseq.design.add_primer(_df[:, "sequence"], 100+i, "rev")
        _df[:, "rev_primer2"] =  fill(100+i, nrow(_df))
    end
    i += 1
end
for _df in gdf[(i - 1)*n_per_group+1:end]
    _df[:, "sequence"] = wgregseq.design.add_primer(_df[:, "sequence"], 100+i, "rev")
    _df[:, "rev_primer2"] =  fill(100+i, nrow(_df))
end


df_final = combine(gdf, :)
df_final.sequence = [ x[1:end-1] for x in df_final.sequence]

if any(length.(df_final.sequence) .!= 231)
    throw(ErrorException("Not all sequences are 231bp!"))
else
    println("Second reverse primer added and last base pair trimmed.")
    println()
end

gdf = groupby(df_final, :genes)

n_per_group = 20
i = 1
while i * n_per_group < length(gdf)
    for _df in gdf[1+(i-1)*n_per_group:i*n_per_group]
        _df[:, "sequence"] = wgregseq.design.add_primer(_df[:, "sequence"], 200+i, "rev")
        _df[:, "rev_primer3"] =  fill(200+i, nrow(_df))
    end
    i += 1
end
for _df in gdf[(i - 1)*n_per_group+1:end]
    _df[:, "sequence"] = wgregseq.design.add_primer(_df[:, "sequence"], 200+i, "rev")
    _df[:, "rev_primer3"] =  fill(200+i, nrow(_df))
end

df_final = combine(gdf, :)
df_final.sequence = [ x[1:end-1] for x in df_final.sequence]

if any(length.(df_final.sequence) .!= 250)
    throw(ErrorException("Not all sequences are 250bp!"))
else
    println("Third reverse primer added and last base pair trimmed.")
    println()
end

##

df_final
## Quality control

