using wgregseq, CSV, DataFrames, SparseArrays, Distributed, SharedArrays
addprocs(6)



@everywhere begin
    function get_data(df)
        seq_vec_0 = vcat([fill(i, df.ct_0[i]) for i in 1:nrow(df)]...)
        mu0 = fill(0, length(seq_vec_0))
    
        seq_vec_1 = vcat([fill(i, df.ct_1[i]) for i in 1:nrow(df)]...)
        mu1 = fill(1, length(seq_vec_1))
    
        mu = vcat(mu0, mu1)
        mu = log10.((df.ct_1 .+ 1) ./ (df.ct_0 .+ 1))
        seq_mat = sparse(vcat([vcat(wgregseq.utils.onehot_encoder.(seq)'...)' for seq in df.promoter]...));
        return seq_mat, mu
    end
    using wgregseq, CSV, DataFrames, SparseArrays, DelimitedFiles, SharedArrays
    
    function run_mcmc(df, i)
        prom = df.name[1]
        # get dataset
        # transform dataset
        seq_mat, mu = get_data(df)
        mat = wgregseq.footprints.run_mcmc(seq_mat, mu, warmup_steps=50000, sample_steps=50000, thin=100)
        writedlm("/Users/tomroeschinger/git/1000_genes_ecoli/code/analysis/20230907_barcode/mcmc_results/$(prom)_$(i)_mcmc.txt", mat)
        return mat
    end
end
    

for i in [1]#[1, 2, 3, 5, 6, 7, 10, 12, 13]
    _df = CSV.read("../../../data/extracted_barcodes/20230907_barcode/$(i)_collapsed.txt", DataFrame)
    gdf = groupby(_df, :name)
    @sync @distributed for _df in gdf
        run_mcmc(_df, i)
    end
end
