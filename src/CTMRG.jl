module CTMRG

using ITensors

include("tensor.jl")
include("tag.jl")
include("diag.jl")
include("1D_ising.jl")
include("2D_ising.jl")
include("spin.jl")
include("CTMRG_prepare.jl")

end # module CTMRG
