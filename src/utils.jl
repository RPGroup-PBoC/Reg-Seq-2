using BioSequences, CSV, DataFrames

# Define custom function for nice imports
Base.parse(::Type{Vector{String}}, x::String) = Vector{String}(filter(x-> x != ", ", split(x, "\""))[2:end-1])
function Base.parse(::Type{Vector{Float64}}, x::String)
    number = split(split(x, "[")[end][1:end-1], ", ")
    number_list = Float64[]
    for num in number
        if num != ""
            push!(number_list, parse(Float64, num))
        else
            return push!(number_list, NaN)
        end
    end
    return number_list

end

function Base.parse(::Type{Tuple{Int64, Tuple{Int64, Int64}}}, x::T) where {T<:String}
    numbers = split(x, ", ")
    ind1 = parse(Int64, split(numbers[1], "(")[2])
    ind2 = parse(Int64, split(numbers[2], "(")[2])
    ind3 = parse(Int64, split(numbers[3], ")")[1])
    return (ind1, (ind2, ind3))
end

Base.parse(::Type{Vector{String}}, x::Missing) = String[]
Base.parse(::Type{Vector{Float64}}, x::Missing) = Float64[]


"""
    onehot_encoder(sequence::BioSequences.LongDNA)

Return One-Hot encoding of DNA sequence.
"""
function onehot_encoder(sequence::BioSequences.LongDNA)

    if occursin(ExactSearchQuery(dna"N"), sequence)
        throw(ArgumentError("Sequence cannot contain 'N'"))
    end
        
    encoder = Dict(DNA_A => 1, DNA_C => 2, DNA_G =>3, DNA_T => 4)
    encoded_seq = zeros(Int64, length(sequence), 4)
    for i in 1:length(sequence)
        encoded_seq[i, encoder[sequence[i]]] = 1
    end
    return encoded_seq
end


"""
    function eval_emat(
        sequence::BioSequences.LongDNA, 
        emat; 
        start_ind_seq::Int64=1, 
        stop_ind_seq::Int64=length(sequence),
        start_ind_emat::Int64=1,
        stop_ind_emat::Int64=size(emat, 1)
    )

Evaluate energy matrix on given sequence.

# Parameters
- sequence: DNA sequence that is to be evaluated. Has to be either a string or BioSequences.LongDNA
- emat: Energy Matrix
- start_ind_seq: First index in the sequence to be evaluated. Default is 1.
- stop_ind_seq: Last index in the sequence to evaluated. Default is the last index of the sequence.
- start_index_seq: First row of energy matrix to be evaluated. Default is 1.
- stop_ind_seq: Last row of energy matrix to be evaluated. Default is the last row of the matrix.

"""
function eval_emat(
        sequence::Union{BioSequences.LongDNA, AbstractString}, 
        emat::Matrix; 
        start_ind_seq::Int64=1, 
        stop_ind_seq::Int64=length(sequence),
        start_ind_emat::Int64=1,
        stop_ind_emat::Int64=size(emat, 1)
    )

    if typeof(sequence) <: AbstractString
        sequence = BioSequences.LongDNA{4}(sequence)
    end
    
    if ismissing(stop_ind_emat)
        stop_ind_emat = size(emat)[1]
    end
    
    if (start_ind_seq < 0) || (start_ind_seq > stop_ind_seq) || (start_ind_seq >length(sequence))
        throw(ArgumentError("Starting index has to be within the sequence. Start Index: $start_ind_seq, Stop Index: $stop_ind_seq"))
    end
    
    if (start_ind_emat < 0) || (start_ind_emat > stop_ind_emat) || (start_ind_emat >size(emat)[1])
        throw(ArgumentError("Starting index has to be within the energy matrix. Start Index: $start_ind_emat, Stop Index: $stop_ind_emat"))
    end
    
    if (stop_ind_seq - start_ind_seq) != (stop_ind_emat-start_ind_emat)
        throw(ArgumentError("Length of sequence not equal to length of energy matrix"))
    end
    
    seq = sequence[start_ind_seq:stop_ind_seq]
    _emat = emat[start_ind_emat:stop_ind_emat, :]
    
    if size(_emat)[1] != length(seq)
        throw(ArgumentError("Energy matrix has to be the same length as sequence"))
    elseif size(_emat)[2] != 4
        throw(ArgumentError("Second dimension has to be of size 4 for energy matrix. Is $(size(emat)[2])"))
    end
            
    onehot_encoder(seq::BioSequences.LongDNA)
    encoded = onehot_encoder(seq)
    
    return sum(encoded .* _emat)
end


"""
    import_twist_order(filename)

Import csv file for twist order. Fixes types in dataframe.
"""
function import_twist_order(filename)
    df = CSV.read(
        filename,
        types=Dict(
            "fwd_primer" => String,
            "rev_primer1" => String,
            "rev_primer2" => String,
            "rev_primer3" => String,
        ),
        DataFrame
    )
    df.fwd_primer = parse.(Tuple{Int64, Tuple{Int64, Int64}}, df.fwd_primer)
    df.rev_primer1 = parse.(Tuple{Int64, Tuple{Int64, Int64}}, df.rev_primer1)
    df.rev_primer2 = parse.(Tuple{Int64, Tuple{Int64, Int64}}, df.rev_primer2)
    df.rev_primer3 = parse.(Tuple{Int64, Tuple{Int64, Int64}}, df.rev_primer3)


    df.genes = parse.(Vector{String}, df.genes)
    df.sequence = [LongDNA{4}(seq) for seq in df.sequence]
    return df
end
