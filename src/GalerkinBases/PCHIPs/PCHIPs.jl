#=
module for piecewise cubic hermite interpolating polynomials
=#
module PCHIPs

export PCHIP

using ..InterpolatingKnots

# Basis functions
# See definition here: https://en.wikipedia.org/wiki/Cubic_Hermite_spline
h00(t) = 2*t.^3 .- 3*t.^2 .+ 1
h10(t) = t.^3 .- 2*t.^2 .+ t
h01(t) = -2*t.^3 .+ 3*t.^2
h11(t) = t.^3 .- t.^2

dh00(t) = 6*t.^2 .- 6*t
dh10(t) = 3*t.^2 .- 4*t .+ 1
dh01(t) = -6*t.^2 .+ 6*t
dh11(t) = 3*t.^2 .- 2*t

d2h00(t) = 12*t .- 6
d2h10(t) = 6*t .- 4
d2h01(t) = -12*t .+ 6
d2h11(t) = 6*t .- 2

# Structure for a pchip spline
mutable struct PCHIP{T}
    knots::Knots
    pvalues::Vector{T}
    mvalues::Vector{T}
end

# Convenience constructor
function PCHIP(x::AbstractVector{Float64}, pvalues::Vector{T}, mvalues::Vector{T}) where T
    return PCHIP(Knots(x), pvalues, mvalues)
end
    
function PCHIP(x::AbstractVector{Float64}, dofs::Vector{T}) where T
    return PCHIP(Knots(x), dofs[1:2:end], dofs[2:2:end])
end

function PCHIP(knots::Knots, dofs::Vector{T}) where T
    return PCHIP(knots, dofs[1:2:end], dofs[2:2:end])
end

function PCHIP(knots::Knots, f::Function, df::Function)
    return PCHIP(knots, f.(knots.x), df.(knots.x))
end

function PCHIP(x::AbstractVector{Float64}, f::Function, df::Function)
    return PCHIP(knots(x), f.(x), df.(x))
end

function PCHIP(knots::Knots, f::PCHIP{T}) where T
    return PCHIP(knots, f.(knots.x), f.(knots.x, 1))
end

function PCHIP(x::AbstractVector{Float64}, f::PCHIP{T}) where T
    return PCHIP(knots(x), f.(x), f.(x, 1))
end

# Evaluate the spline at x in interval k
function (f::PCHIP{T})(k::Int64, x::Real)::T where T
    dx = f.knots.dx[k]
    t = (x - f.knots.x[k])/dx
    return T(h00(t)*f.pvalues[k] + h10(t)*dx*f.mvalues[k] + h01(t)*f.pvalues[k+1] + h11(t)*dx*f.mvalues[k+1] )
end

# Evaluate the derivative of the spline at x in interval k
function (f::PCHIP{T})(k::Int64, x::Real, deriv::Val{1})::T where T
    dx = f.knots.dx[k]
    t = (x - f.knots.x[k])/dx
    return T(dh00(t)*f.pvalues[k]/dx + dh10(t)*f.mvalues[k] + dh01(t)*f.pvalues[k+1]/dx + dh11(t)*f.mvalues[k+1] )
end

# Evaluate the 2nd derivative of the spline at x in interval k
function (f::PCHIP{T})(k::Int64, x::Real, deriv::Val{2})::T where T
    dx = f.knots.dx[k]
    t = (x - f.knots.x[k])/dx
    return T(d2h00(t)*f.pvalues[k]/dx/dx + d2h10(t)*f.mvalues[k]/dx + d2h01(t)*f.pvalues[k+1]/dx/dx + d2h11(t)*f.mvalues[k+1]/dx )
end

# Deal with unknown derivatives?
function (f::PCHIP{T})(k::Int64, x::Real, deriv::Val{N})::T where {T,N}
    return zero(T)
end

# Evaluate the spline at x
function (f::PCHIP{T})(x::Real)::T where T # Return type will be T, i.e., the type of pvalue/mvalue
    k = get_interval(f.knots, x)    
    return f(k, x)
end

function (f::PCHIP{T})(x::Real, deriv::Int64)::T where T # Return type will be T, i.e., the type of pvalue/mvalue
    # Find the interval
    k = get_interval(f.knots, x)
    if k == -1
        return zero(T)
    else
        # Now evaluate in the interval
        dx = f.knots.dx[k]
        t = (x - f.knots.x[k])/dx
        if deriv == 0 # don't calculate a derivative
            return f(k, x)
        else
            return f(k, x, Val(deriv))
        end
    end
end

function Base.vec(f::PCHIP{T}) where T
    ans = Vector{T}(undef, 2*f.knots.N)
    ans[1:2:end] .= f.pvalues
    ans[2:2:end] .= f.mvalues
    return ans
end

include("Matrices.jl")

end # module
