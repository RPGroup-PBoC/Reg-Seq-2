module wgregseq

include("enzyme_list.jl")

module design
include("design.jl")
end

module plotting_style
include("plotting_style.jl")
end

module promoter_finder
include("promoter_finder.jl")
end

end # module
