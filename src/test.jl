module m

using Main.GalerkinBases
using Main.SolveEilenberger

coarseknots = Knots(range(0, stop = 20, length = 11);)
knots = refine_knots(coarseknots, 8)
basis = PCHIPBasis(knots)
fineknots = refine_knots(knots, 16)
x = knots.x
xfine = fineknots.x

# Create ΔΨ's
function zeros_one(n, i)
    out = zeros(n)
    out[i] = 1
    return out
end

n = 2*length(coarseknots)
Δdof = trues(n)
Aydof = trues(n)
Aydof[n-1] = 0
PCHIPzero = PCHIP(knots, zeros(2*length(knots)))
i = 1

δΨs = [[ΨDOF(basis, PCHIP(knots, PCHIP(coarseknots, zeros_one(n,i))), PCHIPzero) for i = 1:n if Δdof[i]];
       [ΨDOF(basis, PCHIPzero, PCHIP(knots, PCHIP(coarseknots, zeros_one(n,i)))) for i = 1:n if Aydof[i]]]

           


# parameters = EilenbergerParameters(1.0, 0.9)
# parameters.Ha = 0.9 * parameters.Hc

# Ψ = ΨDOF(basis, parameters.ΔT)
# solve!(Ψ, parameters)

#=
using LinearAlgebra
n,i,j = 0, 1, 1
Φ = ΦDOF(basis)
ω, ny, nz, ρ0, ρ1, ρ2, ρ3 = ωnρ(parameters, n, i, j)
MΔ_mat = Main.SolveEilenberger.MΔ(Ψ)
MAy_mat = Main.SolveEilenberger.MAy(Ψ)
Main.SolveEilenberger.solve_eilenberger!(Φ, MΔ_mat, MAy_mat, ω, ny, nz, parameters.ΔT)
f = Main.SolveEilenberger.f(Φ)
fb = Main.SolveEilenberger.fb(Φ)
g = Main.SolveEilenberger.g(Φ)

Φf = ΦDOF(basis)
Ψ.dof[1] += 1e-4
MΔ_mat = Main.SolveEilenberger.MΔ(Ψ)
MAy_mat = Main.SolveEilenberger.MAy(Ψ)
Main.SolveEilenberger.solve_eilenberger!(Φf, MΔ_mat, MAy_mat, ω, ny, nz, parameters.ΔT)
Ψ.dof[1] -= 1e-4
f_f = Main.SolveEilenberger.f(Φf)
fb_f = Main.SolveEilenberger.fb(Φf)
g_f = Main.SolveEilenberger.g(Φf)

df_fd(x) = (f_f.(x) - f.(x)) ./1e-4
dfb_fd(x) = (fb_f.(x) - fb.(x)) ./1e-4
dg_fd(x) = (g_f.(x) - g.(x)) ./1e-4

A = Main.SolveEilenberger.A(Φ, MΔ_mat, MAy_mat, ω, ny, nz, parameters.ΔT)
MδΔ_mat = Main.SolveEilenberger.MΔ(δΨ)
MδAy_mat = Main.SolveEilenberger.MAy(δΨ)
δΦ = ΦDOF(basis)
Main.SolveEilenberger.solve_eilenberger_sensitivity!(δΦ, Φ, A, MδΔ_mat, MδAy_mat, ω, ny, nz, parameters.ΔT)

df = Main.SolveEilenberger.f(δΦ)
dfb = Main.SolveEilenberger.fb(δΦ)
dg = Main.SolveEilenberger.g(δΦ)

@show norm(df_fd.(fineknots.x) - df.(fineknots.x), Inf)
@show norm(dfb_fd.(fineknots.x) - dfb.(fineknots.x), Inf)
@show norm(dg_fd.(fineknots.x) - dg.(fineknots.x), Inf)
#=
norm(df_fd.(fineknots.x) - df.(fineknots.x), Inf) = 2.1096591434863615e-12
norm(dfb_fd.(fineknots.x) - dfb.(fineknots.x), Inf) = 2.815151266096748e-12
norm(dg_fd.(fineknots.x) - dg.(fineknots.x), Inf) = 4.03148928473061e-12
=#
=#
end # module
