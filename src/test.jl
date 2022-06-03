module m

using Main.GalerkinBases
using Main.SolveEilenberger

x = range(0, stop = 20, length = 51);
basis = PCHIPBasis(x)

parameters = EilenbergerParameters(1.0, 0.9)

Ψ = ΨDOF(basis, parameters.ΔT)
parameters.Ha = 0.9 * parameters.Hc
solve!(Ψ, parameters)

end # module
