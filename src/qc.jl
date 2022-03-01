using BioSequences, DataFrames, CairoMakie

import ..enzyme_list
using ..wgregseq: design.import_primer

dir = @__DIR__
home_dir = joinpath(split(dir, "/")[1:end-1])

function check_dataframe(df; print_results=true, site_start=1, site_end=nothing)
    gdf = groupby(df, :promoter)
    max_list, min_list = Float64[], Float64[]
    output=""
    for _df in gdf
        promoter = _df.promoter |> unique
        
        if print_results
            println("Promoter: $promoter")
            println("-------------------------")
        else
            output *= "Promoter: $promoter\n"
            output *= "-------------------------\n"
        end

        mut_rate = mutation_coverage(_df, site_start, site_end)
        min_rate, max_rate = minimum(mut_rate), maximum(mut_rate)
        if print_results
            println("Minimum mutation rate: $(min_rate)")
            println("Maximum mutation rate: $(max_rate)")
            println()
            
        else
            output *= "Minimum mutation rate: $(min_rate)\n"
            output *= "Maximum mutation rate: $(max_rate)\n"
            output *= "\n"

        end
        push!(min_list, min_rate)
        push!(max_list, max_rate)
        if check_cut_sites(_df)
            if print_results
                println("Cut sites are correct.")
                println()
            else
                output *= "Cut sites are correct.\n"
                output *= "\n" 
            end
        else
            if print_results
                println("Cut sites NOT are correct.")
                println()
            else
                output *= "Cut sites NOT are correct.\n"
                output *= "\n" 
            end

        end

        if check_primers(_df)
            if print_results
                println("Primers are correct.")
                println()
            else
                output *= "Primers are correct.\n"
                output *= "\n" 
            end
        else
            if print_results
                println("Primes NOT are correct.")
                println()
            else
                output *= "Primes NOT are correct.\n"
                output *= "\n" 
            end

        end

        if print_results
            println()
            println()
        else
            output *= "\n"
            output *= "\n"
        end
    end

    if ~print_results
        open("/$home_dir/qc_output.txt", "w") do file
            write(file, output)
        end
    end
    check_primers_off_target(df)
    promoter = df.promoter |> unique
    fig = Figure(resolution=(15*length(promoter), 800))
    ax = Axis(fig[1, 1])
    
    scatter!(ax, 1:length(promoter), max_list, label="Maximum")
    scatter!(ax, 1:length(promoter), min_list, label="Minimum")
    lines!(ax, [1, length(promoter)], [0.1, 0.1], color="grey", linestyle=:dash)
    ax.ylabel = "Mutation Rate"
    ax.xticklabelrotation = pi/4
    ax.xticks = 1:length(promoter)
    ax.xtickformat = x -> string.(promoter)
    axislegend()
    save("/$home_dir/figures/min_max_mutation_rates.pdf", fig)
    return fig
end





function mutation_coverage(df, site_start=1, site_end=nothing)
    
    if isnothing(site_end)
        site_end = length(df.sequence[1])
    end 

    seq_list = [seq[site_start:site_end] for seq in df.sequence]
    mat = PFM(seq_list)
    mat = mat ./ length(seq_list)
    mut_rate = 1 .- [maximum(A) for A in eachcol(mat)]
    return mut_rate
end


function check_cut_sites(df)
    gdf = groupby(df, ["upstream_re_site", "downstream_re_site"])
    sites_correct = true
    for _df in gdf
        enz1 = unique(_df.upstream_re_site)[1]
        enz2 = unique(_df.downstream_re_site)[1]

        site1 = enzyme_list[enzyme_list.enzyme .== enz1, "site"][1]
        site2 = enzyme_list[enzyme_list.enzyme .== enz2, "site"][1]
        sites_correct *= ~any(map(seq -> seq[21:26] != LongDNA{4}(site1),  _df.sequence))
        sites_correct *= ~any(map(seq -> seq[187:192] != LongDNA{4}(site2),  _df.sequence))
    end
    return sites_correct
end


function check_primers(df)
    primer_cols = filter(x -> occursin("primer", x), names(df))
    gdf = groupby(df, primer_cols)
    primers_correct = true
    for _df in gdf
        for col in primer_cols
            direction = split(col, "_")[1] |> string
            ind, positions  = unique(_df[!, col])[1]

            primer_seq = import_primer(ind, direction)[1:(positions[2]-positions[1]+1)]
            primers_correct *= ~any(map(seq -> seq[positions[1]:positions[2]] != LongDNA{4}(primer_seq),  _df.sequence))
        end
    end
    return primers_correct
end



function check_primers_off_target(df)
    primer_cols = filter(x -> occursin("rev", x), names(df))
    for col in primer_cols
        primers = unique(df[!, col])
        for primer in primers
            primer_len = primer[2][2] - primer[2][1] + 1
            println(primer)
            primer_seq = import_primer(primer[1], "rev")[1:primer_len]
            query = ApproximateSearchQuery(primer_seq)
            inds0 = [~isnothing(findfirst(query, 0, seq)) for seq in df.sequence]
            x0 = (inds0 .!= map(x -> x == primer, df[!, col]))

            x1 = [findfirst(query, 1, seq) for seq in df.sequence]
            x1 = x1[.~ isnothing.(x1)]
            x1 = x1[length.(collect.(x1)) .== primer_len]

            x2 = [findfirst(query, 2, seq) for seq in df.sequence]
            x2 = x2[.~ isnothing.(x2)]
            x2 = x2[length.(collect.(x2)) .== primer_len]

            if sum(x0) > 0
                println("$(sum(x0)) off target binding for primer $(primer[1]) with no mismatches.")
            elseif length(x1) > 0
                println("$(length(x1)) off target binding for primer $(primer[1]) with one mismatch.")
            elseif length(x2) > 0
                println("$(length(x2)) off target binding for primer $(primer[1]) with two mismatches.")
            else
                println("No off target binding for primer $(primer[1]).")
            end
            println()
        end
    end
end

