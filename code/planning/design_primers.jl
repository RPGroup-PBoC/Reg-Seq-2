using wgregseq, DataFrames, CSV, BioSequences, CairoMakie, Statistics, StatsBase


# Set path
dir = @__DIR__
home_dir = joinpath(split(dir, "/")[1:end-2])


##
filename = "2022-02-15_twist_order.csv"

df = CSV.read(
    "/$home_dir/data/twist_orders/$filename",
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

##
primer_df_IDT = DataFrame()
primer_df = DataFrame()

overhang = lowercase(wgregseq.design.import_primer(300, "rev"))

ind0 = 142 
i = 0

enzyme_df = wgregseq.enzyme_list

SbfI_site = enzyme_df[enzyme_df.enzyme .== "SbfI", "site"][1]
SalI_site = enzyme_df[enzyme_df.enzyme .== "SalI", "site"][1]

barcode_IDT =  "(N:25252525)" * join((fill("(N)", 19))) 
barcode =  join((fill("N", 20))) 

for primer in unique(df.fwd_primer)
    primer_ind = primer[1]
    primer_seq = wgregseq.design.import_primer(primer_ind, "fwd")
    name = "SC$(ind0+i)_$(primer_ind)_fwd"
    append!(primer_df, DataFrame(name=name, sequence=primer_seq))
    i += 1
end

for primer in unique(df.rev_primer1)
    primer_ind = primer[1]
    primer_seq = overhang * reverse_complement(wgregseq.design.import_primer(primer_ind, "rev"))
    name = "SC$(ind0+i)_$(primer_ind)_rev"
    append!(primer_df, DataFrame(name=name, sequence=primer_seq))
    i += 1
    #=
    bar_primer_seq = string(reverse_complement(LongDNA{4}(SalI_site))) * 
                        barcode_ * 
                        string(reverse_complement(LongDNA{4}(SbfI_site))) *
                        string(primer_seq)
    bar_primer_name = "SC$(ind0+i)_$(primer_ind)_SbfI_BC_SalI_rev"
    append!(primer_df, DataFrame(name=bar_primer_name, sequence=bar_primer_seq))
    i += 1
    println(bar_primer_seq)
    =#
end

for primer in unique(df.rev_primer2)
    primer_ind = primer[1]
    primer_seq = overhang * reverse_complement(wgregseq.design.import_primer(primer_ind, "rev")[1:19])
    name = "SC$(ind0+i)_$(primer_ind)_rev"
    append!(primer_df, DataFrame(name=name, sequence=primer_seq))
    i += 1
    #=
    bar_primer_seq = string(reverse_complement(wgregseq.design.import_primer(primer_ind, "rev")))
    bar_primer_seq = reverse_complement(LongDNA{4}(bar_primer_seq * SbfI_site * barcode * SalI_site))
    bar_primer_name = "SC$(ind0+i)_$(primer_ind)_SbfI_BC_SalI_rev"
    append!(primer_df, DataFrame(name=bar_primer_name, sequence=string(primer_seq)))
    i += 1
    println(bar_primer_seq)
    =#
end

for primer in unique(df.rev_primer3)
    primer_ind = primer[1]
    primer_seq = overhang * reverse_complement(wgregseq.design.import_primer(primer_ind, "rev")[1:19])
    name = "SC$(ind0+i)_$(primer_ind)_rev"
    append!(primer_df, DataFrame(name=name, sequence=primer_seq))
    i += 1
    #=
    bar_primer_seq = string(reverse_complement(wgregseq.design.import_primer(primer_ind, "rev")))
    bar_primer_seq = reverse_complement(LongDNA{4}(bar_primer_seq * SbfI_site * barcode * SalI_site))
    bar_primer_name = "SC$(ind0+i)_$(primer_ind)_SbfI_BC_SalI_rev"
    append!(primer_df, DataFrame(name=bar_primer_name, sequence=string(primer_seq)))
    i += 1
    println(bar_primer_seq)
    =#
end
spacer = "CGATAAAC"
barcoding_primer_IDT = spacer * string(reverse_complement(LongDNA{4}(SalI_site))) * 
                        barcode_IDT * 
                        string(reverse_complement(LongDNA{4}(SbfI_site))) * 
                        string(overhang)

barcoding_primer = LongDNA{4}(spacer) * reverse_complement(LongDNA{4}(SalI_site)) * 
                        LongDNA{4}(barcode) * 
                        reverse_complement(LongDNA{4}(SbfI_site)) * 
                        overhang
append!(primer_df, DataFrame(name="SC$(ind0+i)_300_SbfI_BC_SalI_rev", sequence=barcoding_primer))
println("Barcoding Primer for IDT: ", barcoding_primer_IDT)
primer_df
CSV.write("/$home_dir/data/twist_orders/primer_table.csv", primer_df)