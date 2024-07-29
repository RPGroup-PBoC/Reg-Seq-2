using BioSequences, CSV, DataFrames
import BioSequences.reverese_complement
using CairoMakie

Base.convert(::Type{String}, x::BioSequences.LongSequence{DNAAlphabet{4}}) = string(x)
Base.convert(::Type{BioSequences.LongDNA}, x::String) = LongDNA{4}(x)

BioSequences.reverse_complement(x::String) = string(reverse_complement(LongDNA{4}(x)))

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
    onehot_encoder(sequence::BioSequences.LongDNA)

Return One-Hot encoding of DNA sequence.
"""
function onehot_encoder(sequence::AbstractString)
    return onehot_encoder(BioSequences.LongDNA{4}(sequence))
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
    #return sum(encoded .* _emat)
end


function eval_emat(
    sequence_list::Union{Vector{T1}, Vector{T2}}, 
    emat::Matrix; 
    start_ind_seq::Int64=1, 
    stop_ind_seq::Int64=length(sequence_list[1]),
    start_ind_emat::Int64=1,
    stop_ind_emat::Int64=size(emat, 1)
) where {T1<:BioSequences.LongDNA, T2<:AbstractString}
    return [eval_emat(
       x, 
        emat; 
        start_ind_seq=start_ind_seq, 
        stop_ind_seq=stop_ind_seq,
        start_ind_emat=start_ind_emat,
        stop_ind_emat=stop_ind_emat
    ) for x in sequence_list]
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


"""
    function joincols(df, column1::Union{Symbol, String}, column2::Union{Symbol, String})

Join two columns of dataframe. If one column has an entry as `missing`, entry from other column
is taken. If both columns have an entry, they need to be identical.
"""
function joincols(df, column1::Union{Symbol, String}, column2::Union{Symbol, String})
    new_col = String[]
    # Iterate through columns
    for (x, y) in zip(df[:, column1], df[:, column2])
        # Check if entries are missing
        miss_arr = [~ismissing(x), ~ismissing(y)]
        if all(miss_arr)
            if x != y
                throw(ErrorException("Columns have non-unique entries"))
            else
                push!(new_col, x)
            end
        elseif any(miss_arr)
            push!(new_col, [x, y][miss_arr][1])
        #elseif all(.~ miss_arr)
           # push!(new_col, missing)
        end
    end
    return new_col
end



num_unique(x) = length(unique(x))


function ecdf!(ax, x)
    lines!(ax, sort(x), 1/length(x):1/length(x):1)
end

function get_mapping_data()
    dir = @__DIR__
    path = joinpath(split(dir, '/')[1:end-1])

    df_map = CSV.read(
        "/$path/data/barcode_maps/20220514_mapping/mapped_barcodes_filtered.csv", 
        DataFrame, 
    )

    # Filter out unnannotad sequences
    df_map = df_map[df_map.name .!= "*", :]

    # Filter out non-unique barcodes
    gdf = groupby(df_map[(df_map.counts .> 2), :], :barcode)
    _df = DataFrame()
    for df in gdf
        if nrow(df) == 1
            append!(_df, df)
        end
    end
    df_map = copy(_df)

    # Get twist order to get wild type sequences
    df_seqs = import_twist_order("/$path/data/twist_orders/2022-02-15_twist_order.csv")
    df_wt = df_seqs[1:1501:119*1501, :];
    insertcols!(df_wt, 4, :promoter_seq => [string(x[27:186]) for x in df_wt.sequence])

    df_wt.promoter_seq |> unique |> length
    df_map = leftjoin(df_map, rename(df_wt[!, [:promoter, :promoter_seq]], :promoter => :name), on="name")
    rename!(df_map, :promoter_seq => :wt_seq);
    return df_map
end


function get_dataset(i, promoter="all";df_map=get_mapping_data())
    dir = @__DIR__
    path = joinpath(split(dir, '/')[1:end-1])

    if  isfile("/$path/data/barcode_counts/20240621_barcode/$(i)_DNA_collapsed.txt")
        df_DNA = CSV.read(
            "/$path/data/barcode_counts/20240621_barcode/$(i)_DNA_collapsed.txt", 
            DataFrame, 
            ignorerepeated=true, 
            delim=" ", 
            header=["ct_0", "barcode"]
        )
        # import RNA
        df_RNA = CSV.read(
            "/$path/data/barcode_counts/20240621_barcode/$(i)_RNA_collapsed.txt", 
            DataFrame, 
            ignorerepeated=true, 
            delim=" ", 
            header=["ct_1", "barcode"]
        )
    elseif isfile("/$path/data/barcode_counts/20231207_barcode/$(i)_DNA_collapsed.txt")
        df_DNA = CSV.read(
            "/$path/data/barcode_counts/20231207_barcode/$(i)_DNA_collapsed.txt", 
            DataFrame, 
            ignorerepeated=true, 
            delim=" ", 
            header=["ct_0", "barcode"]
        )
        # import RNA
        df_RNA = CSV.read(
            "/$path/data/barcode_counts/20231207_barcode/$(i)_RNA_collapsed.txt", 
            DataFrame, 
            ignorerepeated=true, 
            delim=" ", 
            header=["ct_1", "barcode"])
    elseif isfile("/$path/data/barcode_counts/20230907_barcode/$(i)_DNA_collapsed.txt")
        df_DNA = CSV.read(
            "/$path/data/barcode_counts/20230907_barcode/$(i)_DNA_collapsed.txt", 
            DataFrame, 
            ignorerepeated=true, 
            delim=" ", 
            header=["ct_0", "barcode"]
        )
        # import RNA
        df_RNA = CSV.read(
            "/$path/data/barcode_counts/20230907_barcode/$(i)_RNA_collapsed.txt", 
            DataFrame, 
            ignorerepeated=true, 
            delim=" ", 
            header=["ct_1", "barcode"])
    end
    
    # merge DNA and RNA reads
    df = outerjoin(df_DNA, df_RNA, on=:barcode)
    
    # replace missing reads with 0
    replace!(df.ct_0, missing => 0)
    replace!(df.ct_1, missing => 0)
    
    # identify promoter sequences
    df = innerjoin(df, df_map, on=:barcode)
    
    # compute total counts
    insertcols!(df, 1, :ct => df.ct_0 .+ df.ct_1)
    insertcols!(df, 1, :relative_counts => (df.ct_1 .+ 1) ./ (df.ct_0 .+ 1))
    
    # Turn sequences into integer
    insertcols!(df, 3, :int_promoter => make_int.(df[:, :promoter]))
    insertcols!(df, 3, :int_wt => make_int.(df[:, :wt_seq]));
    if promoter != "all"
        df = df[df.name .== promoter, :]
    end
    return df
end


"""
    function get_reps(gc, promoter="all")

Import reads for growth condition `gc` for all replicates.
"""
function get_reps(gc, promoter="all"; df_map=get_mapping_data())

    dir = @__DIR__
    path = joinpath(split(dir, '/')[1:end-1])

    gcs = unique([x[1] for x in split.(vcat(
                                        readdir("/$path/data/barcode_counts/20230907_barcode/"), 
                                        readdir("/$path/data/barcode_counts/20231207_barcode/"), 
                                        readdir("/$path/data/barcode_counts/20240621_barcode/")), "_")])
    gcs = gcs[map(x -> x == gc, [split(y, '-')[1] for y in gcs])]
    df = DataFrame()
    reps = [split(x, '-')[2] for x in gcs]
    for (gc, rep) in zip(gcs, reps)
        _df = get_dataset(gc, promoter; df_map=df_map)
        insertcols!(_df, 4, :replicate => rep)
        append!(df, _df)
    end
    return df
end
