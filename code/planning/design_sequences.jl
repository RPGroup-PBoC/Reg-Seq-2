using wgregseq, FASTX, DataFrames, CSV, BioSequences, CairoMakie, Statistics, StatsBase, Dates

# Set path
dir = @__DIR__
home_dir = joinpath(split(dir, "/")[1:end-2])

# Import genome
re = open(FASTA.Reader, "/$home_dir/data/ecocyc/mg1655_genome.fasta")
wt_sequence = [sequence(record) for record in re][1]


##
# Import gene list to generate sequences for
gene_table = CSV.read("/$home_dir/data/100_genes.csv", DataFrames.DataFrame, comment="#")

# Give IDs to groups for adding primers later
group_dict = Dict(filter(x -> occursin("Antibiotic/toxin", x), gene_table.group) .=> 1)
group_dict["Gold Standard"] = 2
group_dict["Heinemann dataset"] = 3
group_dict["uncharacterized protein"] = 4
group_dict["YmfT_modulon"] = 5
group_dict["YgeV_modulon"] = 6
groups = zeros(nrow(gene_table))
for i in 1:nrow(gene_table)
    if gene_table[i, "group"] in keys(group_dict)
        groups[i] = group_dict[gene_table[i, "group"]]
    else
        groups[i] = 7
    end
end

# Add IDs to table
insertcols!(gene_table, 4, :group_ID => groups)

# Print number of genes per group
println(combine(groupby(gene_table, "group"), nrow))

## Import promoters
# Import promoter list and infer types that can be infered automatically
promoter_list = CSV.read(
    "/$home_dir/data/promoter_list_processed.csv", 
    DataFrames.DataFrame, 
    types=Dict(
        "promoter"=>String,
        "tss"=>Float64,
        "direction"=>String
    )
)

operons_without_promoters = CSV.read(
    "/$home_dir/data/operons_without_promoters.csv", 
    DataFrames.DataFrame, 
    types=Dict(
        "direction"=>String
    )
)



# Replace columns by nicer types
promoter_list.genes = parse.(Vector{String}, promoter_list.genes)
promoter_list.gene_position = parse.(Vector{Float64}, promoter_list.gene_position)

operons_without_promoters.genes = parse.(Vector{String}, operons_without_promoters.genes)
operons_without_promoters.gene_position = parse.(Vector{Float64}, operons_without_promoters.gene_position)

## Some genes may have the wrong synomym
all_gene_list = CSV.read(
    "/$home_dir/data/all_genes_table.csv", 
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

# Replace gene names if they are synonyms
for i in 1:nrow(gene_table)
    name = String(gene_table[i, "name"])
    if name ∉ all_gene_list.gene
        syn = all_gene_list[map(x -> name in x, all_gene_list.synonyms), :gene][1]
        gene_table[i, "name"] = syn
    end
end
gene_group_ID_dict = Dict(gene_table.name .=> gene_table.group_ID)
println("Data imported.")
println()
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

# Create temporary DataFrame
temp_df = DataFrame()

# Go through each gene
for gene_list in df.genes |> unique
    # Find promoters for gene
    sub_df = df[map(x -> x == gene_list, df.genes), :]
    # If more than one promoter, find strongest
    if nrow(sub_df) > 1
        append!(temp_df, DataFrame(wgregseq.design.find_best_promoter(sub_df, wt_sequence)))
    else
        append!(temp_df, sub_df)
    end
end

df = temp_df

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
    "/$home_dir/data/urtecho_2020/tss_operon_regulation.txt", 
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
            #println("No operons for genes $(genes) in Urtecho dataset.")
            #println()
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

# Add found promoters to list
append!(df, df_tss_urtecho[:, ["genes", "tss", "direction", "gene_position", "promoter"]])

# Remove genes with identified promoters from list of genes without promoters
df_no_prom = df_no_prom[Not(delete_index_list), :]
println("Took TSS from Urtecho data set.")
println()

## Now we tackle the genes for which we have not found a TSS yet
# We look for the site with the strongest predicted promoter within
# 500 bases upstream of the start of the coding region
p = wgregseq.promoter_finder.Promoter_Calculator()
tss_list = Int64[]
name_list = String[]

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

    push!(name_list, join(df_no_prom[i, "genes"], "_") * "_predicted")
end

# Add start sites to genes 
insertcols!(df_no_prom, 2, :tss =>tss_list)
# Add information that they are predicted
insertcols!(df_no_prom, 5, :promoter =>name_list)
# Add promoters to list
append!(df, df_no_prom)

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
    global df_sequences = vcat(df_sequences, _df)
end

if any(length.(df_sequences.sequence) .!= 160)
    throw(ErrorException("Not all sequences are 160bp!"))
else
    println("Mutated sequences created!")
    println()
end

## Check sequences for restriction sites 
gdf = groupby(deepcopy(df_sequences), :genes)
df_stack = DataFrame()


enzymes = ["SalI", "SacI", "NheI", "XbaI", "SpeI", "XhoI", "EcoRI", "ApaI", "ScaI", "NcoI", "MluI", "EcoRV", "BbsI", "BamHI", "AgeI"]


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

dict_enz = Dict{Any, Any}(enzymes .=> sum.(eachcol(df_stack[!, enzymes])))
dict_enz["gene"] = [String31["all"]]
dict_enz["promoter"] = "all"
append!(df_stack, DataFrame(dict_enz))
println("")
println(df_stack)

# Set these enzymes to disable query
enz1 = ""
enz2 = ""

println("Upstream restriction enzyme (default is SalI by hitting `enter`):")
while enz1 ∉ wgregseq.enzyme_list.enzyme
    global enz1 = readline()
    if enz1 == ""
        global enz1 = "SalI"
    end
    if enz1 ∉ wgregseq.enzyme_list.enzyme
        println("$enz1 not in list of enzymes")
    end
end
println("Downstream restriction enzyme (default is SacI by hitting `enter`):")
while enz2 ∉ wgregseq.enzyme_list.enzyme
    global enz2 = readline()
    if enz2 == ""
        global enz2 = "SacI"
    end
    if enz2 ∉ wgregseq.enzyme_list.enzyme
        println("$enz2 not in list of enzymes")
    end
end

## Add restriction sites
df_sequences.sequence = wgregseq.design.add_re_sites.(df_sequences.sequence, enz1, enz2)
if any(length.(df_sequences.sequence) .!= 172)
    println(length.(df_sequences.sequence) |> unique)
    throw(ErrorException("Not all sequences are 172 after adding primers bp!"))
end


## Adding primers
# First we check that the primer does not contain the restriction site. If that happens
# we try to take the next primer

primer = wgregseq.design.check_primers_re_sites(enz1, enz2, 100, "both")
insertcols!(df_sequences, 4, :fwd_primer => fill(primer, nrow(df_sequences)))
insertcols!(df_sequences, 5, :rev_primer1 => fill(primer, nrow(df_sequences)))
df_sequences.sequence = wgregseq.design.add_primer(df_sequences.sequence, primer)

# Confirm that the primers are correct
fwd_primer = wgregseq.design.import_primer(primer, "fwd")
rev_primer = wgregseq.design.import_primer(primer, "rev")
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

# Go through groups of n_per_group genes and add primer
i = 1
while i * n_per_group < length(gdf)
    global primer = wgregseq.design.check_primers_re_sites(enz1, enz2, 100 + i, "rev") 
    for _df in gdf[1+(i-1)*n_per_group:i*n_per_group]
        _df[:, "sequence"] = wgregseq.design.add_primer(_df[:, "sequence"], primer, "rev")
        _df[:, "rev_primer2"] =  fill(primer, nrow(_df))
    end
    global i += 1
end
primer = wgregseq.design.check_primers_re_sites(enz1, enz2, primer+1, "rev")
for _df in gdf[(i - 1)*n_per_group+1:end]
    _df[:, "sequence"] = wgregseq.design.add_primer(_df[:, "sequence"], primer, "rev")
    _df[:, "rev_primer2"] =  fill(primer, nrow(_df))
end

# Combine sequences to single table
df_final = combine(gdf, :)
df_final.sequence = [ x[1:end-1] for x in df_final.sequence]

# Confirm sequence length
if any(length.(df_final.sequence) .!= 231)
    throw(ErrorException("Not all sequences are 231bp!"))
else
    println("Second reverse primer added and last base pair trimmed.")
    println()
end

# Group sequences again
gdf = groupby(df_final, :genes)

# Add reverse primer for each group ID that was determined in the beginning
for _df in gdf
    genes = _df[1, "genes"]
    gene_ids = []
    for gene in genes
        if gene in keys(gene_group_ID_dict)
            group = gene_group_ID_dict[gene]
            push!(gene_ids, group)
        end
    end
    gene_ids = gene_ids |> unique

    if length(gene_ids) > 1
        throw(ErrorException("More than one group for genes: $genes."))
    else
        ID = gene_ids[1]
    end
    local primer = wgregseq.design.check_primers_re_sites(enz1, enz2, Int(200 + ID), "rev")
    _df[:, "sequence"] = wgregseq.design.add_primer(_df[:, "sequence"], primer, "rev")
    _df[:, "rev_primer3"] =  fill(primer, nrow(_df))
end

# Combine groups to single table
df_final = combine(gdf, :)
df_final.sequence = [ x[1:end-1] for x in df_final.sequence]

# Confirm Length
if any(length.(df_final.sequence) .!= 250)
    throw(ErrorException("Not all sequences are 250bp!"))
else
    println("Third reverse primer added and last base pair trimmed.")
    println()
end

##
# Save results
filename = string(Dates.today()) * "_sequence_list.csv"
filename_twist = string(Dates.today()) * "_sequence_list_twist.csv"
CSV.write("/$home_dir/data/$filename", df_final)
CSV.write("/$home_dir/data/$filename_twist", df_final[!, "sequence"])
println("Sequence list saved in `/$home_dir/data/$filename`")
println("Total number of sequences: $(nrow(df_final))")

##
println()
println("Sequences per group.")
println(combine(groupby(df_final, :rev_primer3), nrow))


##