mutable struct ΨDOF{T} <: DOF
    basis::GalerkinBasis{T}
    dof::Vector{Float64}
    Aydof::BitVector
end

Δhat(Ψ::ΨDOF) = MethodError(Δhat, (Ψ,)) |> throw
Ayhat(Ψ::ΨDOF) = MethodError(Ayhat, (Ψ,)) |> throw

Δ(Ψ::ΨDOF{T}) where T = T(Ψ.basis.knots, Δhat(Ψ))
Ay(Ψ::ΨDOF{T}) where T = T(Ψ.basis.knots, Ayhat(Ψ))
MΔ(Ψ::ΨDOF) = sum(  Δhat(Ψ) .* basis(Ψ).M)
MAy(Ψ::ΨDOF) = sum(  Ayhat(Ψ) .* basis(Ψ).M)
    
δΨs(coarseknots::Knots, basis::GalerkinBasis{T}) where T = MethodError(δΨs, (coarseknots, basis)) |> throw

################################################################################
# PCHIP Implementation
################################################################################

function ΨDOF(basis::PCHIPBasis, ΔT)
    Ndof = NDOF(basis)
    dof = zeros(Float64, 2*Ndof - 1) 
    N = 2*length(basis.knots)
    dof[1:2:N-1] .= ΔT
    Aydof = trues(Ndof)
    Aydof[end-1] = false # = 0
    return ΨDOF{PCHIP}(basis, dof, Aydof)
end

function ΨDOF(basis::PCHIPBasis, Δ::PCHIP, Ay::PCHIP)
    Ndof = NDOF(basis)
    Aydof = trues(Ndof)
    Aydof[end-1] = false # = 0
    dof = [vec(Δ); vec(Ay)[Aydof]]
    return SolveEilenberger.ΨDOF{PCHIP}(basis, dof, Aydof)    
end

Δhat(Ψ::ΨDOF{PCHIP}) = Ψ.dof[1:NDOF(basis(Ψ))]

function Ayhat(Ψ::ΨDOF{PCHIP})
    out = Vector{Float64}(undef, NDOF(basis(Ψ)))
    out[Ψ.Aydof] = Ψ.dof[NDOF(basis(Ψ))+1:end]
    out[end-1] = 0
    return out
end

# convenience function to create a vector of zeros with a one in spot i
function zeros_one(n, i)
    out = zeros(n)
    out[i] = 1
    return out
end

function δΨs(coarseknots::Knots, basis::PCHIPBasis)
    knots = basis.knots
    n = 2*length(coarseknots)
    Δdof = trues(n)
    Δdof[2] = false # Slope at z = 0
    Δdof[end-1] = false # Value at z = ∞
    Aydof = trues(n)
    Aydof[2] = false # Applied Field
    Aydof[end-1] = false # Value zt z = ∞
    PCHIPzero = PCHIP(knots, zeros(2*length(knots)))

    return [[ΨDOF(basis, PCHIP(knots, PCHIP(coarseknots, zeros_one(n,i))), PCHIPzero) for i = 1:n if Δdof[i]];
            [ΨDOF(basis, PCHIPzero, PCHIP(knots, PCHIP(coarseknots, zeros_one(n,i)))) for i = 1:n if Aydof[i]]]
end
