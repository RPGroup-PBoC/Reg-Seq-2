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

module quality_control
include("qc.jl")
end

module utils
include("utils.jl")
export parse
end


end # module
