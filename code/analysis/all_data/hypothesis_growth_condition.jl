using wgregseq, CairoMakie, Statistics, CSV, DataFrames, StatsBase, Random, KernelDensity

wgregseq.viz.default_makie!()

if length(ARGS) == 0
    throw(ErrorException("No Growth Condition provided."))
end

gc = ARGS[1]


"""
    shuffle_counts!(df)

Shuffle counts in dataframe.
"""

function shuffle_counts(df)
    inds = randperm(nrow(df))
    df.ct_0 = df.ct_0[inds]
    df.ct_1 = df.ct_1[inds]
    df.relative_counts = df.relative_counts[inds]
    df.ct = df.ct[inds]
    return df
end


dir = @__DIR__
path = joinpath(split(dir, '/'))

df_gcs = CSV.read("/$path/growth_conditions.csv", DataFrame)

gc_names = Dict(string.(collect(1:40)) .=> [replace(df_gcs[x, :Condition], " " => "-") for x in collect(1:40)])

df_map = wgregseq.utils.get_mapping_data()


function shuffle_reps(df::AbstractDataFrame, gc_name::AbstractString; d=1, shuffles=1)
    footprints = []
    prom = df.name[1]
    # get dataset
    for gdf in groupby(df, :replicate)  
        if d == 0
            push!(footprints, wgregseq.footprints.mutual_information_mutation(gdf, vec=true))
        else
            x = wgregseq.footprints.mutual_information_mutation(gdf, vec=true)
            push!(footprints, [mean(x[i-d:i+d]) for i in 1+d:160-d])
        end
    end
    # rep1
    data_rep1 = footprints[1] ./ sum(footprints[1])
    #rep2
    data_rep2 = footprints[2] ./ sum(footprints[2])
    # sum of replicates (each replicate normalized)
    data_1d = (footprints[1] ./ sum(footprints[1])) .+ (footprints[2] ./ sum(footprints[2]))
    # difference in replicates over sum of replicates (each replicate normalized)
    data_1d_2 = ((footprints[1] ./ sum(footprints[1])) .- (footprints[2] ./ sum(footprints[2]))) ./ ((footprints[1] ./ sum(footprints[1])) .+ (footprints[2] ./ sum(footprints[2])))
    # 
    data_2d = [[footprints[1][i] / sum(footprints[1]), footprints[2][i] / sum(footprints[2])] for i in 1:length(footprints[1])]
    
    
    shuffle_1d = []
    shuffle_1d_2 = []
    shuffle_2d = []
    shuffle_rep1 = Vector{Float64}[]
    shuffle_rep2 = []

    shuffle_df = deepcopy(df)
    
    # Shuffle data and recompute footprints
    for i in 1:shuffles
        for gdf in groupby(shuffle_df, :replicate)  
            if d == 0
                inds = randperm(nrow(gdf))
                push!(footprints, wgregseq.footprints.mutual_information_mutation_vec(int_wt=gdf.int_wt, int_promoter=gdf.int_promoter, ct=gdf.ct[inds], ct_0=gdf.ct_0[inds], ct_1=gdf.ct_1[inds]))
            else
                inds = randperm(nrow(gdf))
                x = wgregseq.footprints.mutual_information_mutation_vec(int_wt=gdf.int_wt, int_promoter=gdf.int_promoter, ct=gdf.ct[inds], ct_0=gdf.ct_0[inds], ct_1=gdf.ct_1[inds])
                push!(footprints, [mean(x[i-d:i+d]) for i in 1+d:160-d])
            end
        end
        # rep 1
        push!(shuffle_rep1, footprints[2i+1])
        # rep 2
        push!(shuffle_rep2, footprints[2i+2])
        # sum of normalized reps
        push!(shuffle_1d, ((footprints[2i+1] ./ sum(footprints[2i+1])) .+ (footprints[2i+2] ./ sum(footprints[2i+2])))...)
        # difference of reps divided by their sum
        push!(shuffle_1d_2, (((footprints[2i+1] ./ sum(footprints[2i+1])) .- (footprints[2i+2] ./ sum(footprints[2i+2]))) ./ ((footprints[2i+1] ./ sum(footprints[2i+1])) .+ (footprints[2i+2] ./ sum(footprints[2i+2]))))...)
        push!(shuffle_2d, (footprints[2i+1] ./ sum(footprints[2i+1]), (footprints[2i+2] ./ sum(footprints[2i+2]))))
    end
    
    med1 = median(hcat(shuffle_rep1...), dims=2)
    med2 = median(hcat(shuffle_rep2...), dims=2)
    
    shuffle_rep1_mat = hcat(shuffle_rep1...)
    quantiles1 = [quantile(shuffle_rep1_mat[i, :], [0.25, 0.75]) for i in 1:160-2d]
    shuffle_rep2_mat = hcat(shuffle_rep2...)
    quantiles2 = [quantile(shuffle_rep2_mat[i, :], [0.25, 0.75]) for i in 1:160-2d]
    
    shuffle1_flat = vcat((shuffle_rep1 ./ vec(sum(hcat(shuffle_rep1...), dims=1)))...)
    shuffle2_flat = vcat((shuffle_rep2 ./ vec(sum(hcat(shuffle_rep2...), dims=1)))...)
    
    rep1_sub = map(x -> max(x, 0), footprints[1] .- vec(med1))
    rep2_sub = map(x -> max(x, 0), footprints[2] .- vec(med2))
    
    # compute pdfs
    X = vcat([hcat(shuffle_2d[i]...) for i in 1:length(shuffle_2d)]...)'
    B = kde((X[1, :], X[2, :]))
    df_pdf = DataFrame(pdf_KDE=[pdf(B, x...) for x in data_2d], position=collect(-115+d:44-d))
    
    df_pdf.pdf_KDE = [max(0, x) for x in df_pdf.pdf_KDE]
    df_pdf.pdf_KDE .+= minimum(filter(x -> x > 0, df_pdf.pdf_KDE))
    
    ecdf = sort(data_1d)
    cdf = map(x -> sum(shuffle_1d .<= x), ecdf) ./ length(shuffle_1d)
    
    ecdf_2 = sort(data_1d_2)
    cdf_2 = map(x -> sum(shuffle_1d_2 .<= x), ecdf_2) ./ length(shuffle_1d_2)
    
    pdf_cutoff = 1
    colors = df_pdf.pdf_KDE .< pdf_cutoff
    
    ############################
    ### init figure and grid ###
    ############################
    
    fig = Figure(size=(2000, 1000))

    ga = fig[1:4, 1:2] = GridLayout()
    gb = fig[1:4, 3] = GridLayout()
    gc = fig[1:4, 4] = GridLayout()
    gd = fig[5, 1:2] = GridLayout()

    # Plot footprints for data and shuffles
    ax1 = Axis(ga[1, 1], title="$prom $gc_name\nOriginal 1", xlabel="position", ylabel="mutual\n information", xticks=-110:5:40, xticklabelsize=7)
    ax2 = Axis(ga[2, 1], title="Original 2", xlabel="position", ylabel="mutual\n information", xticks=-110:5:40, xticklabelsize=7)
    ax3 = Axis(ga[3, 1], title="Background subtracted Rep 1", xlabel="position", xticks=-110:5:40, xticklabelsize=7)
    ax4 = Axis(ga[4, 1], title="Background subtracted Rep 2", xlabel="position", xticks=-110:5:40, xticklabelsize=7)
    
    barplot!(ax1, -115+d:44-d, footprints[1], color="orange")
    lines!(ax1, -115+d:44-d, vec(med1), color="black")
    lines!(ax1, -115+d:44-d, hcat(quantiles1...)[1, :], color="black", linestyle=:dash)
    lines!(ax1, -115+d:44-d, hcat(quantiles1...)[2, :], color="black", linestyle=:dash)
    
    barplot!(ax2, -115+d:44-d, footprints[2], color="orange")
    lines!(ax2, -115+d:44-d, vec(med2), color="black")
    lines!(ax2, -115+d:44-d, hcat(quantiles2...)[1, :], color="black", linestyle=:dash)
    lines!(ax2, -115+d:44-d, hcat(quantiles2...)[2, :], color="black", linestyle=:dash)
    
    barplot!(ax3, -115+d:44-d, rep1_sub, color="gray")
    barplot!(ax4, -115+d:44-d, rep2_sub, color="gray")
    
    # plot scatter of rep1 vs rep2 and sum of reps vs differenc of reps
    ax11 = Axis(gb[1, 1], xlabel="replicate 1", ylabel="replicate 2", title="Normalized mutual\n information")
    ax22 = Axis(gb[2, 1], xlabel="rep 1 + rep 2", ylabel="(rep 1 - rep 2)/(rep 1 + rep 2)")
    
    # only plot subset of scatters to not overfill plot
    for i in 1:min(shuffles, 10)
        scatter!(ax11, footprints[2i + 1] ./ sum(footprints[2i + 1]), footprints[2i + 2] ./ sum(footprints[2i + 2]), color="gray", label="Shuffle")
        scatter!(ax22, footprints[2i + 1] ./ sum(footprints[2i + 1]) .+ footprints[2i + 2] ./ sum(footprints[2i + 2]), (footprints[2i + 1] ./ sum(footprints[2i + 1]) .- footprints[2i + 2] ./ sum(footprints[2i + 2])) ./ (footprints[2i + 1] ./ sum(footprints[2i + 1]) .+ footprints[2i + 2] ./ sum(footprints[2i + 2])), color="gray", label="Shuffle")
    end
    
    # plot data on top of scatter
    scatter!(ax11, footprints[1] ./ sum(footprints[1]), footprints[2] ./ sum(footprints[2]), color="orange", label="Original")
    scatter!(ax22, footprints[1] ./ sum(footprints[1]) .+ footprints[2] ./ sum(footprints[2]), (footprints[1] ./ sum(footprints[1]) .- footprints[2] ./ sum(footprints[2])) ./ (footprints[1] ./ sum(footprints[1]) .+ footprints[2] ./ sum(footprints[2])), label="Original", color="orange")
    
    
    # compute Kolmogorov-Smirnov test and plot ecdf
    KST = maximum(abs.(collect(1/length(ecdf):1/length(ecdf):1) .- cdf))
    ax = Axis(gc[1, 1], xlabel="rep1 + rep2", ylabel="ECDF")
    lines!(ax, ecdf, 1/length(ecdf):1/length(ecdf):1, label="Data", color="orange")
    lines!(ax, ecdf, cdf, label="Shuffles", color="gray")
    ax.title = "Kolmogorov-Smirnov test: $(round(KST, digits=4))"
    axislegend(ax, position=:rb)
    
    # compute Kolmogorov-Smirnov test and plot ecdf
    KST = maximum(abs.(collect(1/length(ecdf_2):1/length(ecdf_2):1) .- cdf_2))
    ax = Axis(gc[2, 1], xlabel="(rep1 - rep2) / (rep1 + rep2)", title="Kolmogorov-Smirnov test\n (based on rep1 + rep2)", ylabel="ECDF")
    lines!(ax, ecdf_2, 1/length(ecdf_2):1/length(ecdf_2):1, label="Data", color="orange")
    lines!(ax, ecdf_2, cdf_2, label="Shuffles", color="gray")
    ax.title = "Kolmogorov-Smirnov test: $(round(KST, digits=4))"
    axislegend(ax, position=:rb)
    
    
    # compute pdfs
    X = vcat([hcat(shuffle_2d[i]...) for i in 1:length(shuffle_2d)]...)'
    B = kde((X[1, :], X[2, :]))
    df_pdf = DataFrame(pdf_KDE=[pdf(B, x...) for x in data_2d], position=collect(-115+d:44-d))

    pdf_shuffle = [pdf(B, X[1, i], X[2, i]) for i in 1:size(X)[2]]
    m_shuff = mean(pdf_shuffle)
    s_shuff = std(pdf_shuffle)

    df_pdf.pdf_KDE = [max(0, x) for x in df_pdf.pdf_KDE]
    df_pdf.pdf_KDE .+= minimum(filter(x -> x > 0, df_pdf.pdf_KDE))
    
    ax = Axis(gd[1, 1], xlabel="position", ylabel="PDF", xticks=-110:5:40, xticklabelsize=7, yscale=log10)
    scatter!(ax, df_pdf.position, df_pdf.pdf_KDE, label="2D KDE")
    lines!(ax, [-115, 44], [m_shuff, m_shuff], label="2D KDE")
    lines!(ax, [-115, 44], [m_shuff, m_shuff] .- s_shuff, linestyle=:dash)
    axislegend(ax)
    return fig
    
end

df = wgregseq.utils.get_reps(gc; df_map=df_map);

for prom in unique(df.name)
    if prom in ["galEp", "ybeDp2"]
        continue
    end

    if ~ispath("/$path/replicate_test/$prom/")
        mkdir("/$path/replicate_test/$prom/")
    end
    fig = shuffle_reps(df[df.name .== prom, :], gc_names(gc); d=1, shuffles=100)
    save("/$path/replicate_test/$prom/$prom-$(gc_names(gc)).pdf", fig)
    println("$prom done.")
end

println("$gc done")
