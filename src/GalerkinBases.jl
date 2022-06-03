module GalerkinBases

import SparseArrays: sparse, spdiagm

export Knots, NDOF, PCHIPBasis, PCHIP

struct GalerkinBasis{T}
    knots
    G # ∫ u u
    D # ∫ u u'
    H # ∫ u' u'
    M # ∫ u u u
end

NDOF(b::GalerkinBasis) = size(b.uu, 1)



joinpath("GalerkinBases", "InterpolatingKnots.jl") |> include
using .InterpolatingKnots

joinpath("GalerkinBases", "PCHIPBasis.jl") |> include

end # module
