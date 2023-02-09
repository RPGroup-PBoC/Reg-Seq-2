using wgregseq, CSV, DataFrames

df_new = CSV.read("/Volumes/rp_lab_ext/og_regseq_data/promoter_bc_key_per_gene/ykgE_identified.csv", DataFrame)
display(first(df_new, 5))

df_old = CSV.read("/Users/tomroeschinger/git/1000_genes_ecoli/data/regseq_elife/ykgE_mapping.txt", DataFrame, ignorerepeated=true, delim=" ")
sort!(df_old, :ct, rev=true)
display(first(df_old, 5))



println("Total counts in new: ", df_new.count |> sum)
println("Total counts in old: ", df_old.ct |> sum)

println()
i = 2
println(df_old[i, :])
println(df_new[df_new.promoter .== df_old.seq[i], :])
println(df_old[df_old.seq .== df_old.seq[i], :])

