joinpath("PCHIPs", "PCHIPs.jl") |> include
using .PCHIPs
joinpath("PCHIPs", "PLIPs.jl") |> include
using .PLIPs

const PCHIPBasis = GalerkinBasis{PCHIP}

function PCHIPBasis(xknots::AbstractVector)
    knots = InterpolatingKnots.Knots(xknots)
    G = PCHIPs.G(knots)
    D = PCHIPs.D(knots)
    H = PCHIPs.H(knots)
    M = PCHIPs.M(knots)
    return PCHIPBasis(knots, G, D, H, M)
end

NDOF(b::PCHIPBasis) = 2*length(b.knots)
