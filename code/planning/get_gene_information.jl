using CSV, DataFrames, wgregseq, CairoMakie

wgregseq.plotting_style.default_makie!()
##
# Open Gene information
dir = @__DIR__
home_dir = joinpath(split(dir, '/')[1:end-2]...)
s = open("/$home_dir/data/ecocyc/genes.dat") do file
    read(file, String)
end

# Drop Comments
gene_list = split(s, "//")[3:end];
# Separate genes
gene_list = [split(x, '\n') for x in gene_list]
# Separate information
gene_list = [[occursin(" - ", x) ? split(x, " - ") : SubString{String}[] for x in gene] for gene in gene_list]
# Drop empty entries
gene_list = [gene[.~isempty.(gene)] for gene in gene_list]

##
# Make DataFrame

ID_list = String[]
name_list = String[]
TU_list = []
direction_list = String[]

for x in gene_list
    name = filter(x -> x[1] == "UNIQUE-ID", x)
    if ~isempty(name)
        push!(ID_list, name[1][2])
    else
        push!(ID_list, "None")
    end
    
    name = filter(x -> x[1] == "COMMON-NAME", x)
    if ~isempty(name)
        push!(name_list, name[1][2])
    else
        push!(name_list, "None")
    end
    
    direction = filter(x -> x[1] == "TRANSCRIPTION-DIRECTION", x)
    if ~isempty(direction)
        push!(direction_list, direction[1][2])
    else
        push!(direction_list, "none")
    end
    
    components = filter(x -> x[1] == "COMPONENT-OF", x)
    tu_list_gene = []
    if ~isempty(components)
        for component in components
            if occursin("TU", component[2])
                push!(tu_list_gene, component[2])
            end
        end
    else
        push!(tu_list_gene, "none")
    end
    push!(TU_list, tu_list_gene)

    
end
df_genes = DataFrame(ID=ID_list, gene=name_list, direction=direction_list, transcription_units=TU_list)

##
# Open Transcription Units
s = open("/$home_dir/data/ecocyc/transunits.dat") do file
    read(file, String)
end

# Drop Comments
tu_list = split(s, "//")[3:end]
# Separate promoters
tu_list = [split(x, '\n') for x in tu_list]
# Separate information
tu_list = [[occursin(" - ", x) ? split(x, " - ") : SubString{String}[] for x in tu] for tu in tu_list]
# Drop empty entries
tu_list = [tu[.~isempty.(tu)] for tu in tu_list]

##
# Write to DataFrame
ID_list = String[]
promoter_list = []


for x in tu_list
    name = filter(x -> x[1] == "UNIQUE-ID", x)
    if ~isempty(name)
        push!(ID_list, name[1][2])
    else
        push!(ID_list, "None")
    end
    
    components = filter(x -> ((x[1] == "COMPONENTS") && (occursin("PM", x[2]))), x)
    if ~isempty(components)
        for component in components
            if occursin("PM", component[2])
                push!(promoter_list, component[2])
            end
        end
    else
        push!(promoter_list, "none")
    end
end

df_tu = DataFrame(TU_ID=ID_list, promoter_ID=promoter_list)

# Combine Genes and Transcription Units
df_tu.genes = [df_genes[map(x -> tu in x["transcription_units"], eachrow(df_genes)), "gene"] for tu in df_tu.TU_ID]

##
# Open Promoter

s = open("/$home_dir/data/ecocyc/promoters.dat") do file
    read(file, String)
end

# Drop Comments
promoter_list = split(s, "//")[3:end]
# Separate promoters
promoter_list = [split(x, '\n') for x in promoter_list]
# Separate information
promoter_list = [[occursin(" - ", x) ? split(x, " - ") : SubString{String}[] for x in promoter] for promoter in promoter_list]
# Drop empty entries
promoter_list = [promoter[.~isempty.(promoter)] for promoter in promoter_list]

##
# Make DataFrame
ID_list = String[]
name_list = String[]
TSS_list = Float64[]
evidence_list = []

for x in promoter_list
    
    name = filter(x -> x[1] == "UNIQUE-ID", x)
    if ~isempty(name)
        push!(ID_list, name[1][2])
    else
        push!(ID_list, "None")
    end
    
    name = filter(x -> x[1] == "COMMON-NAME", x)
    if ~isempty(name)
        push!(name_list, name[1][2])
    else
        push!(name_list, "None")
    end
    
    TSS = filter(x -> x[1] == "ABSOLUTE-PLUS-1-POS", x)
    if ~isempty(TSS)
        push!(TSS_list, parse(Float64, TSS[1][2]))
    else
        push!(TSS_list, NaN)
    end
    
    cits = filter(x -> x[1] == "CITATIONS", x)
    if ~isempty(cits)
        cit_list = []
        for cit in cits
            if occursin("EV-EXP", cit[2])
                push!(cit_list, "EXP")
            elseif occursin("EV-COMP", cit[2])
                push!(cit_list, "COMP")
            end
            if length(cit_list) == 0 
                push!(cit_list, "none")
            end
        end
        push!(evidence_list, cit_list)
    else
        push!(evidence_list, ["none"])
    end
    
end
df_tss = DataFrame(promoter_ID=ID_list, promoter=name_list, tss=TSS_list, evidence=evidence_list)

##
# Join DataFrames
df_joint = outerjoin(df_tu, df_tss, on = :promoter_ID) |> unique
CSV.write("/$home_dir/data/promoter_list.csv", df_joint)
println("First 20 rows of list:")
println(first(df_joint, 20))

##
# Gold Standard Genes

gold_standards = [
    "rspA",
    "araA",
    "araB",
    "znuC",
    "znuB",
    "xylA",
    "xylF",
    "dicC",
    "relE",
    "relB",
    "ftsK",
    "znuA",
    "lacI",
    "marR",
    "dgoR",
    "ompR",
    "dicA",
    "araC"
]

# Split DataFrame into promoters with TUs and without
df_missing_gene = df_joint[ismissing.(df_joint.genes),:]
df_comp = df_joint[.~ismissing.(df_joint.genes),:]

df_gs = vcat([df_comp[map(x -> gene in x, df_comp.genes), :] for gene in gold_standards]...) |> unique
println()
println("Gold Standard Genes:")
println(df_gs)

CSV.write("/$home_dir/data/gold_standards.csv", df_gs)