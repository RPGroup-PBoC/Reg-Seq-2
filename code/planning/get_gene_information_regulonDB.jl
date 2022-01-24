using CSV, DataFrames, wgregseq, CairoMakie

##

# Open Gene information
dir = @__DIR__
home_dir = joinpath(split(dir, '/')[1:end-2]...)
# Import gene list to generate sequences for
gene_table = CSV.read("/$home_dir/data/100_genes.csv", DataFrame)
##
# Import promoter data
df_DB_prom = CSV.read(
    "/$home_dir/data/regulonDB/promoter.txt", 
    DataFrame, 
    comment="#", 
    header=[
        "PROMOTER_ID",
        "PROMOTER_NAME",
        "PROMOTER_STRAND",
        "POS_1",
        "SIGMA_FACTOR",
        "BASAL_TRANS_VAL",
        "EQUILIBRIUM_CONST",
        "KINETIC_CONST",
        "STRENGTH_SEQ",
        "PROMOTER_SEQUENCE",
        "KEY_ID_ORG",
        "PROMOTER_NOTE",
        "PROMOTER_INTERNAL_COMMENT",
    ]
)

# Remove import artifacts
df_DB_prom = df_DB_prom[.~ ismissing.(df_DB_prom.PROMOTER_NAME), :]

# Change direction to be compatible with other database
direction_dict = Dict("forward" => "+", "reverse" => "-")
df_DB_prom.PROMOTER_STRAND = [direction_dict[x] for x in df_DB_prom.PROMOTER_STRAND]

# Fix type
df_DB_prom.PROMOTER_NAME = convert(Vector{String}, df_DB_prom.PROMOTER_NAME)
df_DB_prom

##
df_DB_genes = CSV.read(
    "/$home_dir/data/regulonDB/gene.txt", 
    DataFrame, 
    comment="#",
    header=[
        "GENE_ID",
        "GENE_NAME",
        "GENE_POSLEFT",
        "GENE_POSRIGHT",
        "GENE_STRAND",
        "GENE_SEQUENCE",
        "GC_CONTENT",
        "CRI_SCORE",
        "GENE_NOTE",
        "GENE_INTERNAL_COMMENT",
        "KEY_ID_ORG",
        "GENE_TYPE",
    ]
)



df_DB_genes = df_DB_genes[.~ ismissing.(df_DB_genes.GENE_NAME), :]
gene_pos = []

# Gene Strand information is missing for some genes, turns out those are split by transposons
df_trans = df_DB_genes[ismissing.(df_DB_genes.GENE_STRAND), :]
df_DB_genes = df_DB_genes[.~ ismissing.(df_DB_genes.GENE_STRAND), :]

for i in 1:nrow(df_DB_genes)
    if df_DB_genes[i, "GENE_STRAND"] == "forward"
        push!(gene_pos, df_DB_genes[i, "GENE_POSLEFT"])
    else
        push!(gene_pos, df_DB_genes[i, "GENE_POSRIGHT"])
    end
end
insertcols!(df_DB_genes, 13, :gene_position => gene_pos)
## Transcription units

df_DB_TU = CSV.read(
    "/$home_dir/data/regulonDB/transcription_unit.txt", 
    DataFrame, 
    comment="#",
    header=[
        "TRANSCRIPTION_UNIT_ID",
        "PROMOTER_ID",
        "TRANSCRIPTION_UNIT_NAME",
        "OPERON_ID",
        "KEY_ID_ORG",
        "TRANSCRIPTION_UNIT_NOTE",
        "TU_INTERNAL_COMMENT",
    ],
    delim='\t'
)

# Remove import artifacts
df_DB_TU = df_DB_TU[.~ ismissing.(df_DB_TU.TRANSCRIPTION_UNIT_ID), :]
# Remove transcription units without promoters
df_DB_TU = df_DB_TU[.~ ismissing.(df_DB_TU.PROMOTER_ID), :]

# Fix type
df_DB_TU.TRANSCRIPTION_UNIT_ID = convert(Vector{String}, df_DB_TU.TRANSCRIPTION_UNIT_ID)
df_DB_TU.PROMOTER_ID = convert(Vector{String}, df_DB_TU.PROMOTER_ID)
df_DB_TU

## Link fof tr

df_DB_link = CSV.read(
    "/$home_dir/data/regulonDB/tu_gene_link.txt", 
    DataFrame, 
    comment="#",
    header=[
        "TRANSCRIPTION_UNIT_ID",
        "GENE_ID",
    ],
    delim='\t'
)
df_DB_link
##
df_DB_tss = DataFrame()
for i in 1:nrow(df_DB_TU)
    # Get TU ID and promoter ID
    TU_ID, promoter_ID= df_DB_TU[i, ["TRANSCRIPTION_UNIT_ID", "PROMOTER_ID"]]
    
    # Find gene IDs of genes in TU
    gene_IDs = df_DB_link[map(x -> x == TU_ID, df_DB_link.TRANSCRIPTION_UNIT_ID), "GENE_ID"]

    # Get gene names and locations
    genes = convert(Vector{String}, df_DB_genes[map(x -> x in gene_IDs, df_DB_genes.GENE_ID), "GENE_NAME"])
    gene_positions = convert(Vector{Int64}, df_DB_genes[map(x -> x in gene_IDs, df_DB_genes.GENE_ID), "gene_position"])

    # get promoter data
    promoter_df = df_DB_prom[df_DB_prom.PROMOTER_ID .== promoter_ID, ["PROMOTER_NAME", "PROMOTER_STRAND", "POS_1",]]
    insertcols!(promoter_df, 2, :genes =>fill(genes, nrow(promoter_df)))
    insertcols!(promoter_df, 3, :gene_position =>fill(gene_positions, nrow(promoter_df)))
    rename!(
        promoter_df,
        Dict(
            "PROMOTER_NAME"=>"promoter",
            "PROMOTER_STRAND"=>"direction",
            "POS_1"=>"tss"
        )
    )
    append!(df_DB_tss, promoter_df)
end
df_DB_tss
CSV.write("/$home_dir/data/promoter_list_regulon_DB.csv", df_DB_tss)
##

df_DB_prom[df_DB_prom.PROMOTER_NAME .== "TSS_1868", :]
filter(x -> "dicA" in x, df_DB_tss.genes)

_d = dropmissing(df_DB_TU, "TRANSCRIPTION_UNIT_NAME")
_d[_d.TRANSCRIPTION_UNIT_NAME .== "dicA", :]

df_DB_TU