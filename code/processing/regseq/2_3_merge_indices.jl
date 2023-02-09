using DataFrames, CSV, BioSequences, wgregseq, CairoMakie, Distances, Statistics

in = string(ARGS[1])
# indeces used in sequences
index_list = [
    "CAAG",
    "AGAG",
    "ATGC",
    "TCTA"
]

# dictionary for growth conditions
gc_dict = Dict(index_list .=> ["LB_DNA", "LB_RNA", "heatshock_DNA", "heatshock_RNA"])
if ~isdir("/Volumes/rp_lab_ext/og_regseq_data/bc_by_gc")
    mkdir("/Volumes/rp_lab_ext/og_regseq_data/bc_by_gc")
end

path = "/Volumes/rp_lab_ext/og_regseq_data/$(in)_bc_by_index/"
filepath = readdir(path; join=true)

for target in values(gc_dict) |> collect
    file="/Volumes/rp_lab_ext/og_regseq_data/bc_by_gc/$target.txt"
    rm(file)
end

for file in filepath
    index = split(split(file, "/")[end], ".")[1]
    if length(index)!=4
        continue
    end
    match = findfirst(x -> hamming(index, x) <= 1, index_list)
    if isnothing(match)
        continue
    end
    target = gc_dict[index_list[match]]
    target_file = "/Volumes/rp_lab_ext/og_regseq_data/bc_by_gc/$target.txt"
    run(pipeline(`cat $file`; stdout="$target_file", append=true))
end


for target in values(gc_dict) |> collect
    file="/Volumes/rp_lab_ext/og_regseq_data/bc_by_gc/$target.txt"
    _df = CSV.read(file, DataFrame, header=["count", "index", "barcode"])
    display(first(_df, 5))
    _df = combine(groupby(_df, ["barcode", "index"]), :count => sum => :count)
    display(first(_df, 5))
    CSV.write(file, _df, delim=" ")
end
