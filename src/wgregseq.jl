module wgregseq

include("enzyme_list.jl")

module promoter_finder
include("promoter_finder.jl")
end

using .promoter_finder
module design
include("design.jl")
end

module plotting_style
include("plotting_style.jl")
end


end # module
