joinpath("PCHIPs", "PCHIPs.jl") |> include
using .PCHIPs
joinpath("PCHIPs", "PLIPs.jl") |> include
using .PLIPs

const PCHIPBasis = GalerkinBasis{PCHIP}

function PCHIPBasis(knots::InterpolatingKnots.Knots)
    G = PCHIPs.G(knots)
    D = PCHIPs.D(knots)
    H = PCHIPs.H(knots)
    M = PCHIPs.M(knots)
    u0 = zeros( 2*length(knots))
    u0[1] = 1
    return PCHIPBasis(knots, G, D, H, M, u0)
end

PCHIPBasis(xknots::AbstractVector) = PCHIPBasis(InterpolatingKnots.Knots(xknots))

NDOF(b::PCHIPBasis) = 2*length(b.knots)
