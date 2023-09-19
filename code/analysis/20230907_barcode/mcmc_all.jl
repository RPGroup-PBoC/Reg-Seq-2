using wgregseq, CSV, DataFrames, CairoMakie, SparseArrays

wgregseq.plotting_style.default_makie!()


df_map = CSV.read(
    "../../../data/barcodes/20220514_mapping/mapped_barcodes.csv", 
    DataFrame, 
);

# Filter out unnannotad sequences
df_map = df_map[df_map.name .!= "*", :]

# Filter out non-unique barcodes
gdf = groupby(df_map[(df_map.map_count .> 2), :], :barcode)
_df = DataFrame()
for df in gdf
    if nrow(df) == 1
        append!(_df, df)
    end
end
df_map = copy(_df);

# Get twist order to get wild type sequences
df_seqs = wgregseq.utils.import_twist_order("../../../data/twist_orders/2022-02-15_twist_order.csv")
df_wt = df_seqs[1:1501:119*1501, :];
insertcols!(df_wt, 4, :promoter_seq => [string(x[27:186]) for x in df_wt.sequence])

df_wt.promoter_seq |> unique |> length
df_map = leftjoin(df_map, rename(df_wt[!, [:promoter, :promoter_seq]], :promoter => :name), on="name")
rename!(df_map, :promoter_seq => :wt_seq)


function get_dataset(i)
    df_DNA = CSV.read(
        "../../../data/extracted_barcodes/20230907_barcode/D$(i)_collapsed.txt", 
        DataFrame, 
        ignorerepeated=true, 
        delim=" ", 
        header=["ct_0", "barcode"]
    )
    # import RNA
    df_RNA = CSV.read(
        "../../../data/extracted_barcodes/20230907_barcode/R$(i)_collapsed.txt", 
        DataFrame, 
        ignorerepeated=true, 
        delim=" ", 
        header=["ct_1", "barcode"]
    )
    
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
    insertcols!(df, 3, :int_promoter => wgregseq.footprints.make_int.(df[:, :promoter]))
    insertcols!(df, 3, :int_wt => wgregseq.footprints.make_int.(df[:, :wt_seq]));
    return df
end

genes = df_map.name |> unique
gcs = [1, 2, 3, 5, 6, 7, 10, 12, 13]

println(genes)
