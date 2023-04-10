using DataFrames, Statistics, StatsBase

"""
    function expression_shift(_df)

Plot promoter footprint and add annotations for annotated sites.
"""
function expression_shift(_df)
    # Create copy
    df = copy(_df)

    # Compute relative (with pseudo counts)
    insertcols!(df, 1, :relative_counts => (df.ct_1 .+ 1) ./ (df.ct_0 .+ 1))
    # Turn sequences into integer
    insertcols!(df, 3, :int_promoter => make_int.(df[:, :promoter]))
    
    freq_mat = frequency_matrix(df)
    # find wild type sequence 
    wt_seq = argmax(sum(freq_mat, dims=1), dims=2) |> vec
    wt_seq = map(x -> x[2], wt_seq)

    function is_mut(x)
        return x .!= wt_seq
    end

    insertcols!(df, 4, :is_mutated => is_mut.(df.int_promoter))
    mean_rel_counts = mean(df.relative_counts)

    ex_shift_arr = zeros(160)
    for (x, seq) in zip(df.relative_counts, df.is_mutated)
        ex_shift_arr[seq] .+= (x - mean_rel_counts) / mean_rel_counts
    end
    ex_shift_arr = ex_shift_arr ./ nrow(df)
    return ex_shift_arr
end


"""
    function mutual_information_bases(_df, nbins=4)

Compute mutual information by binning relative counts and using all bases as variables.
"""
function mutual_information_bases(_df, nbins=4)
    # Create copy
    df = copy(_df)

    # Compute relative (with pseudo counts)
    insertcols!(df, 1, :relative_counts => (df.ct_1 .+ 1) ./ (df.ct_0 .+ 1))
    #df = df[df.relative_counts .> 0.05, :]
    f = fit(Histogram, log.(df.relative_counts))#, nbins=nbins)
    #return f
    bins = f.edges[1] |> collect
    nbins = length(bins) - 1

    insertcols!(df, 1, :bin => map(x -> findfirst(y-> log(x) < y, bins) - 1, df.relative_counts))
    l = length(df.promoter[1])
    # bins, bases
    p = zeros(l, nbins, 4)

    # bins
    for j in 1:l
        counts = countmap([[df.bin[i], df.promoter[i][j]] for i in 1:nrow(df)])
        for key in keys(counts) |> collect
            p[j, key[1], DNA_dict[key[2]]] = counts[key]
        end
    end
    p = p ./ nrow(df)
    mut_information = [sum([clog(p[i, j, k], sum(p[i, :, k]), sum(p[i, j, :])) for j in 1:nbins for k in 1:4]) for i in 1:l]
    return mut_information
end


"""
    function mutual_information_mutation(_df, nbins=4)

Compute mutual information using RNA and DNA counts as well as mutated base indentity.
"""
function mutual_information_mutation(_df)
    # Create copy
    df = copy(_df)

    # mu, m
    l = length(df.promoter[1])
    p = zeros(l, 2, 2)

    # Turn sequences into integer
    insertcols!(df, 3, :int_promoter => make_int.(df[:, :promoter]))
    freq_mat = frequency_matrix(df)

    # find wild type sequence 
    wt_seq = argmax(sum(freq_mat, dims=1), dims=2) |> vec
    wt_seq = map(x -> x[2], wt_seq)
    
    function is_mut(x)
        return x .!= wt_seq
    end

    insertcols!(df, 4, :is_mutated => is_mut.(df.int_promoter))

    for (j, seq) in enumerate(df.is_mutated)
        for i in 1:l
            p[i, 1, seq[i]+1] += df.ct_0[j]
            p[i, 2, seq[i]+1] += df.ct_1[j]
        end
    end
    p ./= sum(df.ct)
    mut_information = [sum([clog(p[i, j, k], sum(p[i, :, k]), sum(p[i, j, :])) for j in 1:2 for k in 1:2]) for i in 1:l]

    return mut_information
end


# transform sequences to integers
DNA_dict = Dict('A' => 1, 'C' => 2, 'G' => 3, 'T' => 4)


"""
    make_int(x, dict=DNA_dict)

Turn DNA sequence into integer sequence.
"""
function make_int(x, dict=DNA_dict)
    return Int64[dict[y] for y in x]
end


"""
    function clog(x, y, z)

Logarithm term if mutual information that can handle zeros.
"""
function clog(x, y, z)
    if x == 0 || y == 0 || z==0
        return 0
    else
        return x * log2(x / (y * z))
    end
end


"""
    function frequency_matrix(df)

Compute frequency matrix for dataframe containing RNA and DNA counts
for sequences in integer format.
"""
function frequency_matrix(df)
    # Create matrix to store frequencies
    freq_mat = zeros(2, 4, 160)

    # Check if sequences exist as integer
    if :int_promoter ∉ names(df)
        int_promoter = make_int.(df.promoter)
    else
        int_promoter = df.int_promoter
    end
    # Sum up total counts of each base
    for (gDNA_counts, cDNA_counts, prom) in zip(df.ct_0, df.ct_1, int_promoter)
        for (i, j) in enumerate(prom)
            freq_mat[1, j, i] += gDNA_counts
            freq_mat[2, j, i] += cDNA_counts
        end
    end
    # Normalize to get frequencies
    freq_mat ./= sum(df.ct)

    return freq_mat
end


function mutual_information_add_model(p::Matrix)
    p1 = sum(p, dims=1)
    p2 = sum(p, dims=2)

    return sum([clog(p[j, i], p1[i], p2[j]) for i in 1:length(p1) for j in 1:length(p2)])
end


function gauge_emat(emat)
    #return emat = (emat .- mean(emat)) ./ std(emat)
    return emat
end