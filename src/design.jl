
using BioSequences, StatsBase, DataFrames, FASTX
import ..enzyme_list
using ..wgregseq: promoter_finder.Promoter_Calculator


"""
    join_seqs(up::String, down::String)

Join two strings.
"""
function join_seqs(up::String, down::String)
    return up * down 
end


"""
    join_seqs(up::BioSequences.LongDNA, down::BioSequences.LongDNA)

Join two DNA sequences.
"""
function join_seqs(up::BioSequences.LongDNA, down::BioSequences.LongDNA)
    return up * down 
end


"""
    add_sites_oligo(seq::BioSequences.LongDNA, up::BioSequences.LongDNA, down::BioSequences.LongDNA)

Add sites upstream and downstream of a DNA sequence.
"""
function add_sites_oligo(seq::BioSequences.LongDNA, up::BioSequences.LongDNA, down::BioSequences.LongDNA)
    return join_seqs(join_seqs(up, seq), down)
end


"""
    add_sites_oligo(seq::BioSequences.LongDNA, up::BioSequences.LongDNA, down::BioSequences.LongDNA)

Add restriction sites for given enzymes upstream and downstream of a DNA sequence.
"""
function add_re_sites(oligo::BioSequences.LongDNA, re1::String, re2::String)
    re1_site = LongDNA{4}(enzyme_list[enzyme_list.enzyme.==re1, "site"][1])
    re2_site = LongDNA{4}(enzyme_list[enzyme_list.enzyme.==re2, "site"][1])
    return add_sites_oligo(oligo, re1_site, re2_site)
end


"""
    import_primer(index::Int, direction::String)

Import a primer with given index from Kosuri collection.
"""
function import_primer(index::Int, direction::String)
    dir = @__DIR__
    path = joinpath(split(dir, '/')[1:end-1])
    if direction == "fwd"
        record = open(FASTA.Reader, "/$path/data/forward_finalprimers.fasta") do r
            collect(r)[index]
        end
        return sequence(record)
    elseif direction == "rev"
        record = open(FASTA.Reader, "/$path/data/reverse_finalprimers.fasta") do r
            collect(r)[index]
        end
        return sequence(record) |> LongDNA{4}|> reverse_complement |> String
    else
        throw(ArgumentError("dir has to be either \"fwd\" or \"rev\""))
    end
end


"""
    add_primer(oligos::Vector{LongSequence{DNAAlphabet{4}}}, index::Int, direction::String="both") 

Add primer for given index to DNA sequence. Can add either forward, reverse or both primers.
"""
function add_primer(oligos::Vector{LongSequence{DNAAlphabet{4}}}, index::Int, direction::String="both")
    if direction == "fwd"
        fwd_primer = import_primer(index, "fwd")
        rev_primer = LongDNA{4}("")
    elseif direction == "rev"
        rev_primer = import_primer(index, "rev")
        fwd_primer = LongDNA{4}("")
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
function find_seq(TSS::Int, strand::String, up::Int, dn::Int, genome::BioSequences.LongDNA)
     
    if strand == "-"
        gene = genome[TSS-dn:TSS+up-1]
    
        outgene = gene |> reverse_complement
        
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



function find_seq(TSS::Int, strand::String, up::Int, dn::Int, genome:: String)
    temp = LongDNA{4}(genome)
    out = find_seq(TSS, strand, up, dn, temp) 
    return string(out[1]), out[2], out[3]
end


function _random_mutation_generator(sequence, rate)::Array{Tuple{Int, Int}, 1}
    num_mutations = Int(floor(length(sequence) * rate))
    positions = sample(1:length(sequence), num_mutations, replace=false)
    mutants = sample(1:3, num_mutations)
    return [(x, y) for (x, y) in zip(positions, mutants)]
end


"""
    random_mutation_generator(sequence, rate, num_mutants)

Create a given number of mutants for a given sequence at a certain rate.
"""
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
    gen_all_single_mutants(sequence; alphabet=[DNA_A, DNA_C, DNA_G, DNA_T])

Generate all single mutations for a given sequence. First sequence in list is the wild type sequence.
"""
function gen_all_single_mutants(sequence; alphabet=[DNA_A, DNA_C, DNA_G, DNA_T])
    mutant_list = [sequence]
    for locus in 1:160
        for mutation in filter(x-> x != sequence[locus], alphabet)
            mut_seq = deepcopy(sequence)
            mut_seq[locus] = mutation
            push!(mutant_list, mut_seq)
        end
    end
    return mutant_list
end

"""
    function mutations_rand(
        sequence::BioSequences.LongDNA{4}, 
        rate::Float64,
        num_mutants::Int;
        site_start=1, 
        site_end=length(sequence)
    )
    

# Parameters
----------
- `sequence` : DNA sequence
- `rate` : rate at which sequences are mutated.
- `num_mutants` : Number of mutated sequences.
- `site_start` : Beginning of the site that is about to be mutated, default 0.
- `site_end` : End of the site that is about to be mutated, default -1 (end).
#Returns
-------
- `mutants` : List of mutant sequences. Each element is a string.
"""
function mutations_rand(
    sequence::BioSequences.LongDNA{4}, 
    rate::Float64,
    num_mutants::Int;
    site_start=1, 
    site_end=length(sequence)
)
        
    mutation_window = sequence[site_start:site_end]
    
    # Create list
    mutants = BioSequences.LongDNA{4}[]
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
function find_restriction_sites(enzyme::String, sequence_list::Vector{LongSequence{DNAAlphabet{4}}})

    ind = findfirst(x -> x == enzyme, enzyme_list.enzyme)
    site = LongDNA{4}(enzyme_list[ind, :site])
    return sum(occursin.((ExactSearchQuery(site),), sequence_list))
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
function find_restriction_sites(enzyme::Vector{String}, sequence_list::Vector{LongSequence{DNAAlphabet{4}}})
    sites = find_restriction_sites.(enzyme, (sequence_list,))
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
wt_sequence : BioSequences.LongDNA
    Sequence of wild type genome.

#Returns
df[ind, :] : Dataframe Row
    Row of DataFrame containing the strongest predicted site.

"""
function find_best_promoter(df::DataFrames.DataFrame, wt_sequence::BioSequences.LongDNA)
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
    check_primers_re_sites(enz1, enz2, primer, direction)

Check if restriction enzymes occur in primer site.
"""
function check_primers_re_sites(enz1, enz2, primer, direction)
    site1 = enzyme_list[enzyme_list.enzyme .== enz1, "site"][1] |> LongDNA{4}
    site2 = enzyme_list[enzyme_list.enzyme .== enz2, "site"][1] |> LongDNA{4}

    clear = false
    if direction == "both"
        while clear == false
            fwd_primer = import_primer(primer, "fwd")
            rev_primer = import_primer(primer, "rev")
            if occursin(ExactSearchQuery(site1), fwd_primer)
                primer += 1
            elseif occursin(ExactSearchQuery(site1), rev_primer)
                primer += 1
            elseif occursin(ExactSearchQuery(site2), fwd_primer)
                primer += 1
            elseif occursin(ExactSearchQuery(site2), rev_primer)
                primer += 1
            else
                clear = true
            end
        end
    elseif direction in ["fwd", "rev"]
        while clear == false
            _primer = import_primer(primer, direction)
            if occursin(ExactSearchQuery(site1), _primer)
                primer += 1
            elseif occursin(ExactSearchQuery(site2), _primer)
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


"""
    function mass_spec_oligo(seq, positions, IDT=false)

Design oligos that can be used for DNA chromotography. A given sequence of random
bases and a PstI restriction site at the 5' of the sense oligo. A tag can be added
that IDT recignizes as a Biotin tag at the 5' end. Returns table of sequences.

# Parameters
------------
- seq: DNA sequence that contains binding site as well as wild type bases at each end
- positions: Bases in the DNA sequence to be mutated for the control
- restriction_site: Enzyme whose cut site is added to the oligo, default `PstI`
- IDT: If `true`, add biotin tag at 5' end that is recognized by IDT. Also changes output to string that can be used as Bulk Input.

# Returns
----------

"""
function mass_spec_oligo(seq::Union{String, LongDNA}, positions::Vector{Int};restriction_site::String="PstI", IDT::Bool=false, name::String=nothing)
    
    if typeof(seq) == LongSequence{DNAAlphabet{4}}
        seq = string(seq)
    end

    # Biotin label
    biotin = "/5Biosg/"

    # Restriction site
    cut_site = enzyme_list[enzyme_list.enzyme.==restriction_site, "site"][1]
    
    # Spacer sequence flanking cut site
    spacer = "ctagct"
    
    # Adapter sequence of cut site and spacers
    adapt = spacer * cut_site * spacer

    # Mutated sense oligo
    mutated_seq = copy(collect(seq))
    for pos in positions
        mutated_seq[pos] = sample(filter(x -> x != mutated_seq[pos], ['A', 'C', 'G', 'T']))
    end

    # Sense sequences
    sense_bind = adapt * seq
    sense_mutated = adapt * join(mutated_seq)
   
    # Antisense sequences
    antisense_bind = reverse_complement(LongDNA{4}(sense_bind)) |> string
    antisense_mutated = reverse_complement(LongDNA{4}(sense_mutated)) |> string

    if IDT
        sense_bind = biotin * sense_bind
        sense_mutated = biotin * sense_mutated
    end

    # Dictionary containing sequences
    if isnothing(name)
        return_dict = Dict(
            "sense binding" => sense_bind, 
            "antisense binding" => antisense_bind,
            "sense mutated" => sense_mutated,
            "antisense mutated" => antisense_mutated
        )
    else
        return_dict = Dict(
            "$(name)_binding" => sense_bind, 
            "$(name)_antisense_binding" => antisense_bind,
            "$(name)_mutated" => sense_mutated,
            "$(name)_antisense_mutated" => antisense_mutated
        )
    end
    
    # Add biotin tags and write output in format for IDT
    if IDT
        if length(antisense_bind) <= 60
            mol = "25nm"
        elseif length(antisense_bind) <= 90
            mol = "100nm"
        elseif length(antisense_bind) <= 100
            mol = "250nm"
        else
            return ErrorException("Oligos have length of $(length(antisense_bind)), but can only by 100bp long.")
        end
        output = ""
        for i in collect(keys(return_dict))
            output *= "$i,$(return_dict[i]),$mol,STD\n"
        end

        return println(output)
    else
        return return_dict
    end
end

# Take array of letters as input
mass_spec_oligo(seq::Vector{Char}, positions::Vector{Int};restriction_site::String="PstI", IDT::Bool=false, name::String=nothing) = mass_spec_oligo(join(seq), positions, IDT=IDT, restriction_site=restriction_site, name=name)