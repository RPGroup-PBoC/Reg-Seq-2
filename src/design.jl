
using BioSequences, StatsBase, DataFrames, FASTX
import ..enzyme_list
using ..wgregseq: promoter_finder.Promoter_Calculator

function join_seqs(up::String, down::String)
    return up * down 
end

function join_seqs(up::BioSequences.LongDNASeq, down::BioSequences.LongDNASeq)
    return up * down 
end


function add_sites_oligo(seq::String, up::String, down::String)
    return 
end

function add_sites_oligo(seq::BioSequences.LongDNASeq, up::BioSequences.LongDNASeq, down::BioSequences.LongDNASeq)
    return join_seqs(join_seqs(up, seq), down)
end


function add_re_sites(oligo::BioSequences.LongDNASeq, re1::String, re2::String)
    re1_site = LongDNASeq(enzyme_list[enzyme_list.enzyme.==re1, "site"][1])
    re2_site = LongDNASeq(enzyme_list[enzyme_list.enzyme.==re2, "site"][1])
    return add_sites_oligo(oligo, re1_site, re2_site)
end


function import_primer(index::Int, direction::String)
    dir = @__DIR__
    path = joinpath(split(dir, '/')[1:end-1])
    if direction == "fwd"
        r = open(FASTA.Reader, "/$path/data/forward_finalprimers.fasta")
        record = collect(r)[index]
        return sequence(record)
    elseif direction == "rev"
        r = open(FASTA.Reader, "/$path/data/reverse_finalprimers.fasta")
        record = collect(r)[index]
        return sequence(record) |> reverse_complement
    else
        throw(ArgumentError("dir has to be either \"fwd\" or \"rev\""))
    end

    
    
end


function add_primer(oligos::Vector{BioSequences.LongDNASeq}, index::Int, direction::String="both")
    if direction == "fwd"
        fwd_primer = import_primer(index, "fwd")
        rev_primer = LongDNASeq("")
    elseif direction == "rev"
        rev_primer = import_primer(index, "rev")
        fwd_primer = LongDNASeq("")
    elseif direction == "both"
        fwd_primer = import_primer(index, "fwd")
        rev_primer = import_primer(index, "rev")
    else
        throw(ArgumentError("direction has to be either \"fwd\", \"rev\" or \"both\"."))
    end
    return add_sites_oligo.(oligos, [fwd_primer], [rev_primer])
end


"""
    find_seq(TSS, strand, up, dn, genome)

Find sequence upstream and downstream of TSS. If gene is reverse 
transcribed, return reverse complement of sequence.

Parameters
----------
TSS : int
    Transcription start site. Refers to position in the genome.
strand : string
    Either '+' for fwd strand or '-' for reverse strand.
up : int
    Number of basepairs upstream of TSS to grab. Relative to TSS according to strand / direction.
dn : int
    Number of basepairs downstream of TSS to grab. Relative to TSS according to strand / direction.
genome : string
    Reference genome as a string 
    
Returns
-------
output : string
    Genetic sequence upstream and downstream of TSS
"""
function find_seq(TSS::Int, strand::String, up::Int, dn::Int, genome::BioSequences.LongDNASeq)
     
    if strand == "-"
        gene = genome[TSS-dn:TSS+up-1]
    
        outgene = LongDNASeq(gene) |> reverse_complement
        
        left_pos = TSS-dn
        right_pos = TSS+up
        
    elseif strand == "+"
        outgene = genome[TSS-up+1:TSS+dn]
        left_pos = TSS-up
        right_pos = TSS+dn
    else
        throw(ArgumentError("`strand` has to be either \"+\" or \"-\""))
    end
    return (outgene, left_pos, right_pos)
end



function _random_mutation_generator(sequence, rate)::Array{Tuple{Int, Int}, 1}
    num_mutations = Int(floor(length(sequence) * rate))
    positions = sample(1:length(sequence), num_mutations, replace=false)
    mutants = sample(1:3, num_mutations)
    return [(x, y) for (x, y) in zip(positions, mutants)]
end


function random_mutation_generator(sequence, rate, num_mutants)
    mutant_list = Vector{Vector{Tuple{Int, Int}}}([])
    for i in 1:num_mutants
        push!(mutant_list, _random_mutation_generator(sequence, rate))
    end
    return mutant_list
end

"""
    function mutate_from_index(sequence, index, alph)

Generate mutated sequence from wild type sequence given a an index and 
an alphabet.

Parameters
----------
- sequence : string
    Wild type seqence
- index : array of tuples
    Each tuple is a position and the mutation at that position.
- alph : array 
    Alphabet to pick mutation from
    
Returns
-------
- str
    Sequence including mutation
"""
function mutate_from_index(sequence, index; alphabet=[DNA_A, DNA_C, DNA_G, DNA_T])
    mut_seq = deepcopy(sequence)
    for (locus, mutation) in index
        mut_seq[locus] = filter(x-> x != mut_seq[locus], alphabet)[mutation]
    end
    return mut_seq
end


"""


Creates single or double mutants.
    
 
#Parameters
----------
- sequence : string
    DNA sequence that is going to be mutated.
- num_mutants : `Int`, default None
    Number of mutant sequences. If None, all possible mutatants are created.
- mut_per_seq : `Int`, default 1
    Number of mutations per sequence.
- site_start : `int``, default 0
    Beginning of the site that is about to be mutated.
- site_end : int, default -1
    End of the site that is about to be mutated.
- alph_type : string, default "DNA"
    Can either be "DNA" for letter sequences, or "Numeric" for integer sequences.
- number_fixed : bool
    If True, the number of mutations is fixed as the rate times length of the sequence.
- keep_wildtype : bool, default False
    If True, adds wild type sequence as first sequence in the list.
    
#Returns
-------
- mutants : list
    List of mutant sequences. Each element is a string.
"""
function mutations_rand(
    sequence::BioSequences.LongDNASeq, 
    rate::Float64,
    num_mutants::Int;
    site_start=1, 
    site_end=length(sequence)
)
        
    mutation_window = sequence[site_start:site_end]
    
    # Create list
    mutants = BioSequences.LongDNASeq[]
    push!(mutants, sequence)
        
    mutant_indices = random_mutation_generator(mutation_window, rate, num_mutants)
    
    for x in mutant_indices
        push!(mutants, sequence[1:site_start-1] * mutate_from_index(mutation_window, x) * sequence[site_end+1:end])
    end
    return mutants
end


"""
    find_restriction_sites(enzyme, sequence_list)

Searches for restriction sites of a specific enzyme in a list of sequences
    
#Parameters
----------
enzyme : string
    Name of the enzyme.
sequence_list: array-type
    list of suquences
    
"""
function find_restriction_sites(enzyme::String, sequence_list::Vector{LongDNASeq})

    ind = findfirst(x -> x == enzyme, enzyme_list.enzyme)
    site = LongDNASeq(enzyme_list[ind, :site])
    return sum(occursin.((site,), sequence_list))
end


"""
    find_restriction_sites(enzyme, sequence_list)

Searches for restriction sites of a specific enzyme in a list of sequences
    
#Parameters
----------
enzyme : string
    Name of the enzyme.
sequence_list: array-type
    list of suquences
    
"""
function find_restriction_sites(enzyme::Vector{String}, sequence_list::Vector{LongDNASeq})
    sites = find_restriction_sites.(enzyme, (sequence_list::Vector{LongDNASeq},))
    return DataFrame(enzyme=enzyme, sites=sites)

end


"""
    function find_best_promoter(df)

Find the promoter with the lowest predicited free energy, using the model from
La Fleur et al, 2022.

#Parameters
-----------
df : DataFrame
    DataFrame containing the transcription start sites for promoters.
wt_sequence : BioSequences.LongDNASeq
    Sequence of wild type genome.

#Returns
df[ind, :] : Dataframe Row
    Row of DataFrame containing the strongest predicted site.

"""
function find_best_promoter(df::DataFrames.DataFrame, wt_sequence::BioSequences.LongDNASeq)
    if "tss" ∉ names(df)
        throw(ArgumentError("DataFrame has no column \"tss\"."))
    end
    if "direction" ∉ names(df)
        throw(ArgumentError("DataFrame has no column \"tss\"."))
    end
    # Find range of possible start sites
    min_tss = df.tss |> minimum
    max_tss = df.tss |> maximum

    # Get sequence around start sites
    sequence = wt_sequence[Int(min_tss)-115:Int(max_tss)+115]
    
    # Run model on sequence
    p = Promoter_Calculator()
    r = p(sequence)
    if df.direction[1] == "-"
        sites = [r["Reverse_Predictions_per_TSS"][x - (Int(min_tss)-115)] for x in df.tss]
    elseif df.direction[1] == "+"
        sites = [r["Forward_Predictions_per_TSS"][x - (Int(min_tss)-115)] for x in df.tss]
    else
        throw(ArgumentError("Direction has to be either \"+\" or \"-\"."))
    end
    # Find strongest site
    ind = argmin([site["dG_total"] for site in sites])
    return df[ind, :]
end


"""
"""
function check_primers_re_sites(enz1, enz2, primer, direction)
    site1 = enzyme_list[enzyme_list.enzyme .== enz1, "site"][1] |> LongDNASeq
    site2 = enzyme_list[enzyme_list.enzyme .== enz2, "site"][1] |> LongDNASeq

    clear = false
    if direction == "both"
        while clear == false
            fwd_primer = import_primer(primer, "fwd")
            rev_primer = import_primer(primer, "rev")
            if occursin(site1, fwd_primer)
                primer += 1
            elseif occursin(site1, rev_primer)
                primer += 1
            elseif occursin(site2, fwd_primer)
                primer += 1
            elseif occursin(site2, rev_primer)
                primer += 1
            else
                clear = true
            end
        end
    elseif direction in ["fwd", "rev"]
        while clear == false
            _primer = import_primer(primer, direction)
            if occursin(site1, _primer)
                primer += 1
            elseif occursin(site2, _primer)
                primer += 1
            else
                clear = true
            end
        end
    else
        throw(ArgumentError("direction has to be in [\"both\", \"fwd\", \"rev\"]."))
    end
    return primer
end