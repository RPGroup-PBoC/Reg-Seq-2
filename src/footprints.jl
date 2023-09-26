using DataFrames, Statistics, StatsBase, LinearAlgebra
import AbstractMCMC
import StatsBase
using Distributions
using Random
import MCMCChains
using KernelDensity
using SparseArrays
using ..wgregseq: utils.onehot_encoder

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

    # Turn sequences into integer
    if "int_promoter" ∉ names(df)
        insertcols!(df, 3, :int_promoter => make_int.(df[:, :promoter]))
    end

    if "int_wt" ∉ names(df)
        insertcols!(df, 3, :int_wt => make_int.(df[:, :wt_seq]))
    end
    if :wt_seq ∉ names(df)
        freq_mat = frequency_matrix(df)[1]
        # find wild type sequence 
        wt_seq = argmax(freq_mat, dims=2) |> vec
        wt_seq = map(x -> x[2], wt_seq)

        function is_mut(x)
            return x .!= wt_seq
        end

        insertcols!(df, 4, :is_mutated => is_mut.(df.int_promoter))
    else
        function is_mut(x, y)
            return x .!= y
        end
        insertcols!(df, 4, :is_mutated => is_mut.(df.int_promoter, df.int_wt))
    end

    # initiate distribution
    p = zeros(l, 2, 2)
    
    # helper function for broadcasting
    s(x) = .~(x)

    p[:, 1, 2] = sum(df.is_mutated .* df.ct_0, dims=1)[1]
    p[:, 1, 1] = sum( s.(df.is_mutated) .* df.ct_0, dims=1)[1]
    p[:, 2, 2] = sum(df.is_mutated .* df.ct_1, dims=1)[1]
    p[:, 2, 1] = sum(s.(df.is_mutated) .* df.ct_1, dims=1)[1]

    # normalize
    p ./= sum(df.ct)

    # compute mutual information
    mut_information = [sum([clog(p[i, j, k], sum(p[i, :, k]), sum(p[i, j, :])) for j in 1:2 for k in 1:2]) for i in 1:l]
    
    return mut_information
end


# transform sequences to integers
DNA_dict = Dict('A' => 1, 'C' => 2, 'G' => 3, 'T' => 4)
DNA_dict_rev = Dict(1 => 'A', 2 => 'C', 3 => 'G',4 => 'T')


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
    #=
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
    =#
    ct_0 = sum(onehot_encoder.(df[!, :promoter]) .* df[!, :ct_0]) / sum(df[!, :ct])
    ct_1 = sum(onehot_encoder.(df[!, :promoter]) .* df[!, :ct_1]) / sum(df[!, :ct])
    return ct_0, ct_1
end

"""
    function mutual_information_add_model(p::Matrix)

Compute mutual information for joint probability distribution of two random
variables.
"""
function mutual_information_add_model(p::Matrix)
    p1 = sum(p, dims=1)
    p2 = sum(p, dims=2)

    return sum([clog(p[j, i], p1[i], p2[j]) for i in 1:length(p1) for j in 1:length(p2)])
end

"""
    function gauge_emat(emat)

Apply gauge to matrix. Each column sums to zero and total matrix norm
is 1.
"""
function gauge_emat(emat)
    emat = emat .- sum(emat, dims=1) / 4
    emat = emat ./ norm(emat)
    return emat
end


"""
    mutable struct MetropolisHastings{T, D, C} <: AbstractMCMC.AbstractSampler

Sampler construct. Contains initial parameter guess, proposal distribution and variance of
proposal distribution.
"""
mutable struct MetropolisHastings{T, D, C} <: AbstractMCMC.AbstractSampler 
    init_θ::T
    proposal::D
    sigma::C
end
MetropolisHastings(init_θ::Matrix{<:Real}) = MetropolisHastings(init_θ, reshape(MvNormal(zero(vec(init_θ)), I), 4, 160), 1.)


"""
    struct DensityModel{F<:Function} <: AbstractMCMC.AbstractModel

Model construct. Contains density function, sequences in One-Hot code, and count indentity array.
"""
struct DensityModel{F<:Function} <: AbstractMCMC.AbstractModel
    ℓπ::F
    mu_arr::AbstractArray
    seq_mat
end


"""
    struct Transition{T, L}

Transition construct. Created at each step, contains parameter and log-probability.
"""
struct Transition{T, L}
    θ::T
    lp::L
end
Transition(model::DensityModel, θ) = Transition(θ, ℓπ(model, θ))


"""
    function AbstractMCMC.step(
        rng::AbstractRNG,
        model::DensityModel,
        spl::MetropolisHastings,
        kwargs...
    )  

Initial step function. Returns transition with initial guess.
"""
function AbstractMCMC.step(
    rng::AbstractRNG,
    model::DensityModel,
    spl::MetropolisHastings,
    kwargs...
)   
    return Transition(model, spl.init_θ), spl
end


"""
    function AbstractMCMC.step(
        rng::AbstractRNG,
        model::DensityModel,
        spl::MetropolisHastings,
        θ_prev_T::Transition,
        ;
        kwargs...
    ) 

Step function. Guesses new step and computes density at new proposal. Compute acceptance
probability and return new Transition if step is acceptence. Return previous transition
otherwise.
"""
function AbstractMCMC.step(
    rng::AbstractRNG,
    model::DensityModel,
    spl::MetropolisHastings,
    θ_prev_T::Transition,
    ;
    kwargs...
)
    # Generate a new proposal.
    θ_T = propose(spl, model, θ_prev_T)
    # Calculate the log acceptance probability.
    α = ℓπ(model, θ_T) - ℓπ(model, θ_prev_T) + q(spl, θ_prev_T, θ_T) - q(spl, θ_T, θ_prev_T)

    # Decide whether to return the previous θ or the new one.
    #if log(rand(rng)) < min(α, 0.0)
    if log(rand(rng)) < min(α, 0.0)
        return θ_T, spl, true
    else
        return θ_prev_T, spl, false
    end
end


# Define a function that makes a basic proposal depending on a univariate
# parameterization or a multivariate parameterization.
propose(spl::MetropolisHastings, model::DensityModel, θ::Matrix{<:Real}) = 
    Transition(model, gauge_emat(θ + rand(spl.proposal)))
propose(spl::MetropolisHastings, model::DensityModel, t::Transition) =
    propose(spl, model, t.θ)

# Calculates the probability `q(θ|θcond)`, using the proposal distribution `spl.proposal`.
q(spl::MetropolisHastings, θ::Matrix{<:Real}, θcond::Matrix{<:Real}) =
    logpdf(spl.proposal, θ - θcond)
q(spl::MetropolisHastings, t1::Transition, t2::Transition) = q(spl, t1.θ, t2.θ)

# Calculate the density of the model given some parameterization.
ℓπ(model::DensityModel, θ) = model.ℓπ(model.seq_mat, model.mu_arr, θ)
ℓπ(model::DensityModel, T::Transition) = model.ℓπ(model.seq_mat, model.mu_arr, T.θ)


"""
    function density(seq_mat, mu::Vector{Float64}, θ::Matrix{Float64})

Compute kerned density estimation for additive model and count indentities. Compute
mutual information given KDE.
"""
function density(seq_mat, mu::Vector{Float64}, θ::Matrix{Float64})
    en = (seq_mat * vec(θ))
    y = kde((en, mu), npoints=(10, 512))
    y.density ./= sum(y.density)
    return mutual_information_add_model(y.density) * 100000
end


"""
    function AbstractMCMC.bundle_samples(
        ℓ::DensityModel, 
        s::MetropolisHastings, 
        N::Integer, 
        ts::Vector{<:Transition},
        chain_type::Type{Any};
        param_names=missing,
        kwargs...
    )

Bundle samples into array given an array of transitions.
"""
function AbstractMCMC.bundle_samples(
    ts::Vector{<:Transition},
    param_names=missing,
    kwargs...
)
    # Turn all the transitions into a vector-of-vectors.
    vals = copy(reduce(hcat,[vcat(vec(t.θ), t.lp) for t in ts])')

    # Check if we received any parameter names.
    if ismissing(param_names)
        param_names = ["Parameter $i" for i in 1:(length(first(vals))-1)]
    end

    # Add the log density field to the parameter names.
    push!(param_names, "lp")

    # Bundle everything up and return a Chains struct.
    return vals, param_names
end


"""
    function adapt_sigma(rate)

Adapt variance of proposal distribution given current acceptance rate. Taken from PyMC3.
"""
function adapt_sigma(rate)
    if rate < 0.001
        return 0.1
    elseif rate < 0.05
        return 0.5
    elseif rate < 0.2 
        return 0.9
    elseif rate > 0.5
        return 1.1
    elseif rate > 0.75
        return 2
    elseif rate > 0.95
        return 10
    else
        return 1
    end
end


"""
    function StatsBase.sample(
        model::AbstractMCMC.AbstractModel,
        sampler::AbstractMCMC.AbstractSampler,
        nsamples::Integer,
        seq_mat,
        rng::Random.AbstractRNG=Random.default_rng(),
        adapt_steps::Integer=1000,
        ;
        kwargs...
    )

Sample function. Takes model and sampler constructs. 
"""
function StatsBase.sample(
    model::AbstractMCMC.AbstractModel,
    sampler::AbstractMCMC.AbstractSampler,
    nsamples::Integer,
    seq_mat,
    ;
    rng::Random.AbstractRNG=Random.default_rng(),
    adapt_steps::Integer=1000,
    thin::Integer=100,
    kwargs...
)
    # Obtain the initial sample and state.
    sample, sampler = AbstractMCMC.step(rng, model, sampler; kwargs...)

    ## Save the sample.
    samples = AbstractMCMC.samples(sample, model, sampler, nsamples; kwargs...)
    samples = AbstractMCMC.save!!(samples, sample, 1, model, sampler, nsamples; kwargs...)
    # Step through the sampler.
    acceptance = 0
    for i in 2:nsamples
        # Obtain the next sample and state.
        sample, sampler, accept = AbstractMCMC.step(rng, model, sampler, sample; kwargs...)
        acceptance += accept
        # Save the sample.
        if i%thin == 0
            samples = AbstractMCMC.save!!(samples, sample, i, model, sampler, nsamples; kwargs...)
            println(density(seq_mat, model.mu_arr, sample.θ))
        end
        if i%1000 ==0
            println("$i of $nsamples done.")
        end
        if i%adapt_steps == 0
            sampler.sigma *= adapt_sigma(acceptance/adapt_steps)
            acceptance = 0
            sampler.proposal = reshape(MvNormal(zeros(640), I * sampler.sigma), 4, 160)
        end
    end

    return AbstractMCMC.bundle_samples(samples; kwargs...)
end


"""
    function run_mcmc(
        seq_mat,
        mu;
        warmup_steps=50000,
        sample_steps=450000,
        density=density,
        )

Run MCMC on dataset.
"""
function run_mcmc(
    seq_mat,
    mu;
    warmup_steps=50000,
    sample_steps=450000,
    density=density,
    thin=100,
    adapt_steps=1000
    )
    if warmup_steps%thin != 0
        throw(ArgumentError("Give warmup steps as multiple of thin."))
    elseif sample_steps%thin != 0
        throw(ArgumentError("Give iter steps as multiple of thin."))
    end

    total_steps = warmup_steps + sample_steps

    # Construct a DensityModel.
    model = DensityModel(density, mu, seq_mat)

    # Set up our sampler with initial parameters.
    spl = MetropolisHastings(randn(4, 160))

    # Run sampler
    chain = StatsBase.sample(model, spl, total_steps, seq_mat, thin=thin, adapt_steps=adapt_steps)

    # Return parameters
    x = reshape(mean(chain[1][Int64(warmup_steps/thin):Int64(total_steps/thin), 1:640], dims=1), 4, 160)
end