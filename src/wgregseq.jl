module wgregseq

include("enzyme_list.jl")

module utils
include("utils.jl")
export parse,num_unique,convert,ecdf!
end
using .utils

module promoter_finder
include("promoter_finder.jl")
end

using .promoter_finder
module design
include("design.jl")
end

module viz
include("viz.jl")
end

module quality_control
include("qc.jl")
end

module footprints
include("footprints.jl")
end



end # module
