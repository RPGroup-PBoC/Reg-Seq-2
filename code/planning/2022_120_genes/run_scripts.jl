println("Processing Ecocyc...")
include("../get_gene_information_ecocyc.jl")
println("Done!\n")
println("Processing RegulonDB...")
include("../get_gene_information_regulonDB.jl")
println("Done!\n")
println("Processing search results...")
include("../database_processing.jl")
println("Done!\n")
println("Designing sequences...")
include("design_sequences.jl")
println("Done!\n")
