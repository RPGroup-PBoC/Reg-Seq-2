using Pandas, CairoMakie, DataFrames, LinearAlgebra, Statistics, wgregseq, Jedi

Jedi.styles.default_makie!()

# Import data. Change path if necessary
df = DataFrames.DataFrame(Pandas.read_excel("data/heinemann_data.xlsx", sheet_name="Table S6", header=[2]))
df[!, "Annotated functional COG group (description)"] |> unique
df.Gene |> unique |> length
##
growth_conditions = ["Glucose",
 "LB",
 "Glycerol + AA",
 "Acetate",
 "Fumarate",
 "Glucosamine",
 "Glycerol",
 "Pyruvate",
 "Chemostat µ=0.5",
 "Chemostat µ=0.35",
 "Chemostat µ=0.20",
 "Chemostat µ=0.12",
 "Stationary phase 1 day",
 "Stationary phase 3 days",
 "Osmotic-stress glucose",
 "42°C glucose",
 "pH6 glucose",
 "Xylose",
 "Mannose",
 "Galactose ",
 "Succinate",
 "Fructose"]

 # Compute mean copy number per growth condition
tot_prot = [sum(df[.~isnan.(df[!, cond]), cond]) for cond in growth_conditions]

# Plot total protein number
fig = Figure(resolution = (600, 400))
ax1 = Axis(
    fig[1, 1], 
    ylabel="Total number of proteins",
    xlabel="Growth Conditions"
)
ind = sortperm(tot_prot, rev=true)
lines!(ax1, 1:22, tot_prot[ind])
scatter!(ax1, 1:22, tot_prot[ind], )

ax1.xticklabelrotation = pi/8
ax1.xticks = 1:length(growth_conditions)
ax1.xtickformat = x -> growth_conditions[ind]
ax1.xminorgridwidth = 0

CairoMakie.save("figures/heinemann_total_protein.pdf", fig)
fig

## Function to plot protein number per gene
"""
    copy_number_plot(genes, normalized=False)  

Plot the copy number of a list of genes in all 
growth conditions.
"""
function copy_number_plot(genes, normalized=false)

    fig = Figure(resolution = (600, 400))
    ax1 = Axis(
        fig[1, 1], 
        ylabel="Total number of proteins",
        xlabel="Growth Conditions"
    )

    # Iterate through genes
    for gene in genes
        # Extract row for gene
        if ~(gene in df[!, "Gene"])
            throw(ArgumentError("Gene $gene is not in the data set."))
        end
        row = df[df[!, "Gene"] .== gene, growth_conditions] |> Matrix |> vec
        mask = map(x -> ~isnan(x), row)

        # Normalize
        if normalized
            row[mask] .= row[mask] ./ tot_prot[mask]
        end
        
        # Plot and make label
        lines!(ax1, collect(1:22)[mask], row[mask], label=gene)
        scatter!(ax1, collect(1:22)[mask], row[mask])
    end


    # Annotate plot
    if normalized
        ax1.ylabel = "Normalized Protein Copy Number"
    else
        ax1.ylabel = "Protein Copy Number"
    end
    
    ax1.xticklabelrotation = pi/8
    ax1.xticks = 1:length(growth_conditions)
    ax1.xtickformat = x -> growth_conditions
    ax1.xminorgridwidth = 0
    axislegend(ax1, position= :lt)

    return fig
end

fig = copy_number_plot(["galK", "galM"], true)
fig

##
function compute_variation(df)
    gene_vals = Matrix(df[!, growth_conditions])
    replace!(gene_vals, NaN=>0)
    return var(gene_vals ./ vec(tot_prot)', dims=2) ./ mean(gene_vals ./ vec(tot_prot)', dims=2)
end

df.variation = vec(compute_variation(df))
sort!(df, :variation, rev=true)
copy_number_plot(df.Gene[1:5])
## Uncharacterized genes
unchar_genes = [
    "yabP",
    "yabQ",
    "yacC",
    "yacH",
    "yadG",
    "yadI",
    "yadE",
    "yadM",
    "yadN",
    "yadS",
    "yaeR",
    "yafT",
    "yaeQ",
    "yafZ",
    "ykfM",
    "yagJ",
    "yagK",
    "yagM",
    "yagL",
    "ykgJ",
    "ykgL",
    "ykgR",
    "ykgE",
    "ykgH",
    "yahC",
    "yahL",
    "yahM",
    "yahN"
]

fig = copy_number_plot(
    filter(x -> x in df.Gene, unchar_genes), true
    )
CairoMakie.save("figures/heineman_unidentified.pdf", fig)
##

df_poor = df[ (df[!, "Annotated functional COG class"] .== "POORLY CHARACTERIZED") .& 
(df[!,"Peptides.used.for.quantitation"] .> 5), :]



#ga = copy_number_plot(df_poor.Gene[1:4])
#copy_number_plot(df_poor.Gene[5:8])
#copy_number_plot(df_poor.Gene[9:12])
#copy_number_plot(df_poor.Gene[13:16])
#copy_number_plot(df_poor.Gene[17:20])


##
ref_genes = [
    "cysG",
    "hcaT",
    "idnT",
    "rssA"
]

fig = copy_number_plot(
    filter(x -> x in df.Gene, ref_genes), true
    )
fig


## Compare two genes individually

function compare_genes(gene, growth_condition1, growth_condition2)
    if growth_condition1 ∉ growth_conditions
        throw(ArgumentError("Growth conditions $growth_condition1 not found!"))
    elseif growth_condition2 ∉ growth_conditions
        throw(ArgumentError("Growth conditions $growth_condition2 not found!"))
    end
    if (length(df[df[!, "Gene"] .== gene, growth_condition1]) > 0) && (length(df[df[!, "Gene"] .== gene, growth_condition2]) > 0)
        protein1 = df[df[!, "Gene"] .== gene, growth_condition1][1]
        protein2 = df[df[!, "Gene"] .== gene, growth_condition2][1]
        return protein1 / protein2 * sum(df[!, growth_condition2]) / sum(df[!, growth_condition1])
    else
        return NaN
    end
end
g1 = "LB"
g2 = "Stationary phase 1 day"
x = filter(x -> x > 0, compare_genes.(df[!, "Gene"], g1, g2))
fig = Figure(resolution=(400, 400))
ax = Axis(fig[1, 1], xscale = log10)
lines!(ax, sort(x), range(0, stop=1, length=length(x)))
ax.xlabel = "Fold Change"
ax.ylabel = "ECDF"
ax.title = "Fold Change of normalized protein copy number \n $g1/$g2"

fig