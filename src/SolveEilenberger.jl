module SolveEilenberger

export ΨDOF, ΦDOF, MΔ, MAy, update!
export EilenbergerParameters, ωnρ
# export solve_eilenberger!
export solve!

using Main.GalerkinBases
using NLsolve

include(joinpath("Eilenberger", "DOFs.jl"))
include(joinpath("Eilenberger", "EilenbergerParameters.jl"))
include(joinpath("Eilenberger", "solve_eilenberger.jl"))
include(joinpath("Eilenberger", "self_consistent_equations.jl"))

end # module
