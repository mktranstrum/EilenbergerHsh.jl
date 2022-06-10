module m

using Main.GalerkinBases
using Main.SolveEilenberger

coarseknots = Knots(range(0, stop = 20, length = 11);)
knots = refine_knots(coarseknots, 8)
basis = PCHIPBasis(knots)
fineknots = refine_knots(knots, 16)
x = knots.x
xfine = fineknots.x

parameters = EilenbergerParameters(1.0, 0.9)
parameters.Ha = 0.9 * parameters.Hc

δΨs = SolveEilenberger.δΨs(coarseknots, basis)

Ψ = ΨDOF(basis, parameters.ΔT)
solve!(Ψ, parameters)


using LinearAlgebra
n,i,j = 0, 1, 1
Φ = ΦDOF(basis)
ω, ny, nz, ρ0, ρ1, ρ2, ρ3 = ωnρ(parameters, n, i, j)
Main.SolveEilenberger.solve_eilenberger!(Φ, MΔ(Ψ), MAy(Ψ), ω, ny, nz, parameters.ΔT)
f = Main.SolveEilenberger.f(Φ)
fb = Main.SolveEilenberger.fb(Φ)
g = Main.SolveEilenberger.g(Φ)

Φf = ΦDOF(basis)
h = sqrt(eps())
for (i, δΨ) = δΨs |> enumerate
    @show "Perturbation $i"
    Ψ.dof += h .* δΨ.dof
    Main.SolveEilenberger.solve_eilenberger!(Φf, MΔ(Ψ), MAy(Ψ), ω, ny, nz, parameters.ΔT)
    Ψ.dof -= h .* δΨ.dof
    f_f = Main.SolveEilenberger.f(Φf)
    fb_f = Main.SolveEilenberger.fb(Φf)
    g_f = Main.SolveEilenberger.g(Φf)

    df_fd(x) = (f_f.(x) - f.(x)) ./h
    dfb_fd(x) = (fb_f.(x) - fb.(x)) ./h
    dg_fd(x) = (g_f.(x) - g.(x)) ./h

    A = Main.SolveEilenberger.A(Φ, MΔ(Ψ), MAy(Ψ), ω, ny, nz, parameters.ΔT)
    δΦ = ΦDOF(basis)
    Main.SolveEilenberger.solve_eilenberger_sensitivity!(δΦ, Φ, A, MΔ(δΨ), MAy(δΨ), ω, ny, nz, parameters.ΔT)

    df = Main.SolveEilenberger.f(δΦ)
    dfb = Main.SolveEilenberger.fb(δΦ)
    dg = Main.SolveEilenberger.g(δΦ)

    @show norm(df_fd.(fineknots.x) - df.(fineknots.x), Inf)
    @show norm(dfb_fd.(fineknots.x) - dfb.(fineknots.x), Inf)
    @show norm(dg_fd.(fineknots.x) - dg.(fineknots.x), Inf)
end

end # module
