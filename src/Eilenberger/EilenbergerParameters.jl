import FastGaussQuadrature: gausslegendre

mutable struct EilenbergerParameters
    κ₀
    t
    T
    ΔT
    Hc
    ρ
    N
    m
    Ha
    k    
    nij
    cosθs
    weights
end

function EilenbergerParameters(κ₀, t, ρ = uniformρ; N = 1, m = 2)
    T = Tc*t
    ΔT = calculateΔT(N, T)
    Hc = calculateHc(κ₀, N, T, ΔT)
    cosθs, weights = gausslegendre(m)
    return EilenbergerParameters(κ₀, t, T, ΔT, Hc, ρ, N, m, 0.0, 0.0, Iterators.product(0:N, 1:div(m,2), 1:m), cosθs, weights)
end

const Tc = 0.567

const inv4π = 1/(4π)

uniformρ(θ, ϕ) = inv4π

function ωnρ(parameters::EilenbergerParameters, n, i, j)
    ω = 2π*parameters.T*(n + 0.5)
    nz = parameters.cosθs[i]
    θ = acos(nz)
    ϕ = j*π/parameters.m
    ny = sin(θ)*sin(ϕ)
    # Weights in Fermi-surface integral
    ρ0 = parameters.ρ(θ, ϕ)*parameters.weights[i]*π/parameters.m # Primary point
    ρ1 = parameters.ρ(π-θ, ϕ+π)*parameters.weights[i]*π/parameters.m # Antipodal point
    ρ2 = parameters.ρ(π-θ, ϕ)*parameters.weights[i]*π/parameters.m # Equator reflection
    ρ3 = parameters.ρ(θ, ϕ+π)*parameters.weights[i]*π/parameters.m # Antipodal of equator reflection
    return ω, ny, nz, ρ0, ρ1, ρ2, ρ3
end

# Functions for calculating zero-field values
import ForwardDiff: derivative

function ΔTresidual(ΔT, N, T)
    ΔT2 = ΔT^2
    ans = 0
    for n = 0:N
        omega = 2*pi*T*(n + 0.5)
        ans += 1/omega - 1/sqrt(omega^2 + ΔT2)
    end
    ans *= 2*pi*T
    ans += log(T/Tc)
    return ans
end

function calculateΔT(N, T; tol = 1e-12)
    guess = 0.8
    val = ΔTresidual(guess, N, T)
    while abs(val) > tol
        dval = derivative(DT->ΔTresidual(DT, N, T), guess)
        guess -= val/dval
        val = ΔTresidual(guess, N, T)
    end
    return guess
end

function calculateHc(κ₀, N, T, ΔT)
    ΔT2 = ΔT^2
    ans = ΔT2*log(T/Tc)
    for n = 0:N
        omega = 2*pi*T*(n + 0.5)
        ans += ΔT2/omega - 2*(ΔT2 + omega^2)/sqrt(ΔT2 + omega^2) + 2*omega
    end
    return sqrt(-3*ans/κ₀)
end
