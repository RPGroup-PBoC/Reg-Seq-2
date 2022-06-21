using CSV, CairoMakie, DataFrames, LinearAlgebra, Statistics, wgregseq
wgregseq.plotting_style.default_makie!()

# Import data. Change path if necessary
#using Pandas
#Pandas.to_csv(Pandas.read_excel("data/heinemann_data.xlsx", sheet_name="Table S6", header=[2]), "data/heinemann_data_cleaned.csv")
df = CSV.read("data/heinemann_data_cleaned.csv", DataFrames.DataFrame)
df = df[.~ismissing.(df.Gene), :]
df[!, "Annotated functional COG group (description)"] |> unique
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
tot_prot = [sum(df[.~ismissing.(df[!, cond]), cond]) for cond in growth_conditions]

# Plot total protein number
fig = Figure(resolution = (600, 400))
ax1 = Axis(
    fig[1, 1], 
    ylabel="Total number of proteins",
    xlabel="Growth Conditions"
)
# Sort growth conditions by total protein count
ind = sortperm(tot_prot, rev=true)
growth_conditions = growth_conditions[ind]

# Plot total protein number
lines!(ax1, 1:22, tot_prot[ind])
scatter!(ax1, 1:22, tot_prot[ind], )

ax1.xticklabelrotation = pi/8
ax1.xticks = 1:length(growth_conditions)
ax1.xtickformat = x -> growth_conditions
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
        mask = map(x -> ~ismissing(x), row)

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
    gene_vals = coalesce.(gene_vals, 0)
    return var(gene_vals ./ vec(tot_prot)', dims=2) ./ mean(gene_vals ./ vec(tot_prot)', dims=2)
end

df.variation = vec(compute_variation(df))
sort!(df, :variation, rev=true)
copy_number_plot(df.Gene[1:5])
unknown_ind = coalesce.(df[!, "Annotated functional COG group (description)"] .== "Function unknown", false)
df[unknown_ind, "Gene"][1:20]
copy_number_plot(df[unknown_ind, "Gene"][1:20], true)
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

## Fold change of gene expression
df = df[.~ismissing.(df.Gene), :]

# Compute fold change
function get_all_fold_changes(df, gene, growth_conditions=growth_conditions)
    total_per_cond = [sum(df[.~ismissing.(df[!, cond]), cond]) for cond in growth_conditions]
    mask =  .~ ismissing.(collect(eachrow(df[df.Gene .== gene, growth_conditions])[1]))
    rel_abundance = collect(eachrow(df[df.Gene .== gene, growth_conditions])[1])[mask] ./ total_per_cond[mask]
    return rel_abundance' ./ rel_abundance, mask
end

# Choose gene and fold change limit
gene = "tolC"
lim = 2

mat, mask = get_all_fold_changes(df, gene)
fig = Figure(resolution=(32*sum(mask), 23*sum(mask)))
ax = Axis(fig)

heatmap!(ax, log.(mat), colorrange=(-log2(lim), log2(lim)))

ax.xticks = 1:sum(mask)
ax.xtickformat = x -> growth_conditions[mask]
ax.xminorgridwidth = 0
ax.xticklabelrotation = pi/8
ax.yticks = 1:sum(mask)
ax.ytickformat = x -> growth_conditions[mask]
ax.yminorgridwidth = 0
ax.yreversed = true
fig[1, 1] = ax
Colorbar(fig[1, 2], limits = (-log2(lim), log2(lim)), colormap = :viridis, label = "Log2 Fold Change")
ax.title = gene
fig
CairoMakie.save("figures/heineman_$(gene)_comp_to_all_fc$lim.pdf", fig)


##

inds_max = findall(x -> x .> lim, mat)
y_unique = unique([x[2] for x in inds_max])
x_unique = unique([x[1] for x in inds_max])
mat_reduced = mat[x_unique, y_unique]


fig = Figure(resolution=(700,500))
ax = Axis(fig, aspect = length(x_unique)/length(y_unique))

heatmap!(ax, log.(mat_reduced), colorrange=(-log2(lim), log2(lim)))

ax.xticks = 1:length(x_unique)
ax.xtickformat = x -> growth_conditions[mask][x_unique]
ax.xminorgridwidth = 0
ax.xticklabelrotation = pi/8
ax.yticks = 1:length(y_unique)
ax.ytickformat = x -> growth_conditions[mask][y_unique]
ax.yminorgridwidth = 0
fig[1, 1] = ax
Colorbar(fig[1, 2], limits = (-log2(lim), log2(lim)), colormap = :viridis, label = "Log2 Fold Change")
ax.title = gene
fig
CairoMakie.save("figures/heineman_$(gene)_comp_to_all_fc$(lim)_reduced.pdf", fig)


## Now we choose subsets of the growth conditions and repeat

sub_growth_conditions = ["Glucose",
 "LB",
 "Glycerol + AA",
 "Chemostat µ=0.5",
 "Chemostat µ=0.35",
 "Chemostat µ=0.20",
 "Chemostat µ=0.12",
 "Stationary phase 1 day",
 "Stationary phase 3 days",
 "Osmotic-stress glucose",
 "42°C glucose",
 "pH6 glucose"]

 carbon_growth_conditions = [
    "Glucose",
    "Acetate",
    "Fumarate",
    "Glucosamine",
    "Glycerol",
    "Pyruvate",
    "Xylose",
    "Mannose",
    "Galactose ",
    "Succinate",
    "Fructose"
    ]

function fold_change_to_any(df, growth_conditions)
    cutoff = 2:10
    rel_identified = Float64[]

    for i in cutoff
        identified = 0

        for gene in df.Gene
            mat = get_all_fold_changes(df, gene, growth_conditions)[1]
            mat = coalesce.(mat, 0)
            identified += any(mat .> i)
        end

        push!(rel_identified, identified / length(df.Gene))
    end
    return rel_identified
end

fig = Figure(resolution=(350, 350))
ax = Axis(fig)

rel_identified_all = fold_change_to_any(df, growth_conditions)
rel_identified_sub = fold_change_to_any(df, sub_growth_conditions)
rel_identified_carbon = fold_change_to_any(df, carbon_growth_conditions)

lines!(ax, cutoff, rel_identified_all, linewidth=2, label="All conditions")
lines!(ax, cutoff, rel_identified_sub, linewidth=2, label="Only Glucose as Carbon Source")
lines!(ax, cutoff, rel_identified_carbon, linewidth=2, label="Only varying Carbon")

axislegend()
ax.xlabel = "Fold Change"
ax.ylabel = "Identified Genes"
ax.yticks = 0:0.2:1
ax.title = "Compared to any \nother growth condition"
ylims!(ax, (-0.02, 1.02))
fig[1, 1] = ax
fig
CairoMakie.save("figures/heineman_identified_genes_comp_to_all.pdf", fig)


## Fold change of gene expression relative to median across conditions
df = df[.~ismissing.(df.Gene), :]
function get_fold_changes(df, gene, growth_conditions=growth_conditions)
    total = [sum(df[.~ismissing.(df[!, cond]), cond]) for cond in growth_conditions]
    mask = .~ ismissing.(collect(eachrow(df[df.Gene .== gene, growth_conditions])[1]))
    rel_abundance = collect(eachrow(df[df.Gene .== gene, growth_conditions])[1])[mask] ./ total[mask]
    return rel_abundance ./ median(rel_abundance), mask
end
gene = "galK"
get_fold_changes(df, gene)
fig = Figure(resolution=(600, 300))
ax = Axis(fig)

y, mask = get_fold_changes(df, gene)
lines!(ax, y, xticklabels=growth_conditions)


ax.xticks = 1:sum(mask)
ax.xtickformat = x -> growth_conditions[mask]
ax.xminorgridwidth = 0
ax.xticklabelrotation = pi/8
fig[1, 1] = ax
ax.title = gene
fig
CairoMakie.save("figures/heineman_$(gene)_rel_to_median.pdf", fig)
##

function fold_change_to_median(df, growth_conditions)
    cutoff = 2:10
    rel_identified = Float64[]
    for i in cutoff
        identified = 0
        for gene in df.Gene
            mat = get_fold_changes(df, gene, growth_conditions)[1]
            mat = coalesce.(mat, 0)
            identified += any((mat .> i) .| (1/i .> mat))
        end
        push!(rel_identified, identified / length(df.Gene))
    end
    return rel_identified
end

fig = Figure(resolution=(350, 350))
ax = Axis(fig)

rel_identified_all = fold_change_to_median(df, growth_conditions)
rel_identified_sub = fold_change_to_median(df, sub_growth_conditions)
rel_identified_carbon = fold_change_to_median(df, carbon_growth_conditions)

lines!(ax, cutoff, rel_identified_all, linewidth=2, label="All conditions")
lines!(ax, cutoff, rel_identified_sub, linewidth=2, label="Only Glucose as Carbon Source")
lines!(ax, cutoff, rel_identified_carbon, linewidth=2, label="Only varying Carbon")
axislegend()

ax.xlabel = "Fold Change"
ax.ylabel = "Identified Genes"
ax.title = "Compared to the median \nacross growth conditions"
ax.yticks = 0:0.2:1
ylims!(ax, (-0.02, 1.02))

fig[1, 1] = ax
fig
CairoMakie.save("figures/heineman_identified_genes_rel_to_median.pdf", fig)

## Find 


