mutable struct ΦDOF{T} <: DOF
    basis::GalerkinBasis{T}
    dof::Vector{Complex{Float64}}
    sdof::BitVector
    ddof::BitVector
    gdof::BitVector
    fL::Float64
    gL::Float64
end

function update!(Φ::ΦDOF, fL, gL)
    Φ.fL = fL
    Φ.gL = gL
    nothing
end
    
fL(Φ::ΦDOF) = Φ.fL
gL(Φ::ΦDOF) = Φ.gL

shat(Φ::ΦDOF{T}) where T = MethodError(shat, (Φ,)) |> throw
dhat(Φ::ΦDOF{T}) where T = MethodError(dhat, (Φ,)) |> throw
ghat(Φ::ΦDOF{T}) where T = MethodError(ghat, (Φ,)) |> throw

fhat(Φ::ΦDOF) = 0.5 * (shat(Φ) + dhat(Φ))
fbhat(Φ::ΦDOF) = 0.5 * (shat(Φ) - dhat(Φ))

f(Φ::ΦDOF{T}) where T = T(basis(Φ).knots, fhat(Φ))
fb(Φ::ΦDOF{T}) where T = T(basis(Φ).knots, fbhat(Φ))
g(Φ::ΦDOF{T}) where T = T(basis(Φ).knots, ghat(Φ))


################################################################################
# PCHIP Implementation
################################################################################

function ΦDOF(basis::PCHIPBasis)
    Ndof = NDOF(basis)
    dof = Vector{Complex{Float64}}(undef, 3*Ndof - 4)
    sdof = BitVector([i != Ndof - 1 for i = 1:Ndof]) # s_indices
    ddof = BitVector([i ∉ [1, Ndof - 1] for i = 1:Ndof ]) # d_indices
    gdof = BitVector([i != Ndof - 1 for i = 1:Ndof]) # g_indices
    return ΦDOF{PCHIP}(basis, dof, sdof, ddof, gdof,0,0)
end

function shat(Φ::ΦDOF{PCHIP})
    Ndof = NDOF(basis(Φ))
    out = Vector{Complex{Float64}}(undef, Ndof)
    out[ Φ.sdof] .= Φ.dof[1:Ndof-1]
    out[ Ndof-1 ] = 2*fL(Φ)
    return out
end

function dhat(Φ::ΦDOF{PCHIP})
    Ndof = NDOF(basis(Φ))
    out = Vector{Complex{Float64}}(undef, Ndof)
    out[ Φ.ddof] .= Φ.dof[Ndof:2*Ndof-3]
    out[ 1] = 0
    out[ Ndof - 1] = 0
    return out
end

function ghat(Φ::ΦDOF{PCHIP})
    Ndof = NDOF(basis(Φ))
    out = Vector{Complex{Float64}}(undef, Ndof)
    out[ Φ.gdof] .= Φ.dof[2*Ndof-2:3*Ndof - 4]
    out[ Ndof - 1] = gL(Φ)
    return out
end

