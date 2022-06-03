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

Δhat(Ψ::ΨDOF{PCHIP}) = Ψ.dof[1:NDOF(basis(Ψ))]

function Ayhat(Ψ::ΨDOF{PCHIP})
    out = Vector{Float64}(undef, NDOF(basis(Ψ)))
    out[Ψ.Aydof] = Ψ.dof[NDOF(basis(Ψ))+1:end]
    out[end-1] = 0
    return out
end
