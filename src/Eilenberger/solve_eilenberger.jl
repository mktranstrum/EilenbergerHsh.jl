import SparseArrays: spzeros

spzeros(mask::BitVector) = spzeros(sum(mask))
spzeros(mask1::BitVector, mask2::BitVector) = spzeros(sum(mask1), sum(mask2))
spzeros(mask1::BitVector, m::Int) = spzeros(sum(mask1), m)

sparse_zeros = spzeros

function fL_gL(ω, ΔT)
    d = sqrt(ω^2 + ΔT^2)
    return ΔT/d, ω/d
end


"""
Solve the Eilenberger equations
"""
function solve_eilenberger!(Φ::ΦDOF{PCHIP}, MΔ, MAy, ω, ny, nz, ΔT)
    update!(Φ, fL_gL(ω, ΔT)...)
    update!(Φ,
            M(Φ, MΔ, MAy, ω, ny, nz, ΔT) \ collect(B(Φ, MΔ, MAy, ω, ny, nz, ΔT))
            )
    nothing
end

function M(Φ::ΦDOF{PCHIP}, MΔ, MAy, ω, ny, nz, ΔT)
    b = basis(Φ)
    sdof = Φ.sdof
    ddof = Φ.ddof
    gdof = Φ.gdof
    return [  (nz*b.D)[sdof, sdof]  (ω*b.G - im*ny*MAy)[sdof, ddof]  spzeros(sdof, gdof);
              (ω*b.G - im*ny*MAy)[ddof, sdof]  (nz*b.D)[ddof, ddof]  (-2*MΔ)[ddof, gdof];
              spzeros(gdof, sdof)  (-MΔ)[gdof, ddof]  (2*nz*b.D)[gdof, gdof] ]
end

function B(Φ::ΦDOF{PCHIP}, MΔ, MAy, ω, ny, nz, ΔT)
    b = basis(Φ)    
    sdof = Φ.sdof
    ddof = Φ.ddof
    gdof = Φ.gdof
    fL = Φ.fL
    gL = Φ.gL
    L = NDOF(b) - 1
    return [-nz*b.D[sdof, L]*2*fL;
           -(ω*b.G - im*ny*MAy)[ddof, L]*2*fL + (MΔ)[ddof, L]*2*gL;
            -2*nz*b.D[gdof, L]*gL]
end
