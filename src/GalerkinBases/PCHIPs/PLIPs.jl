#=
module for piecewise linear interpolating polynomials
=#
module PLIPs

export PLIP

using ..InterpolatingKnots

# Basis functions
h0(t) = 1 - t
h1(t) = t

dh0(t) = -1
dh1(t) = 1

# Structure for a plip spline
mutable struct PLIP{T<:Real}
    knots::Knots
    pvalues::Vector{T}
end

# Convenience constructor
function PLIP(knots::AbstractVector{Float64}, pvalues::Vector{T}) where T <: Real
    return PLIP(Knots(knots), pvalues)
end
    
# Evaluate the spline at x
# TODO: This can be made more efficient for evaluating arrays
function (f::PLIP{T})(x::Real, deriv::Int64 = 0)::T where T <: Real # Return type will be T, i.e., the type of pvalue/mvalue
    # Find the interval
    k = get_interval(f.knots, x)
    if k == -1
        return zero(T)
    else
        # Now evaluate in the interval
        dx = f.knots.dx[k]
        t = (x - f.knots.x[k])/dx
        if deriv == 0 # don't calculate a derivative
            return T(h0(t)*f.pvalues[k] + h1(t)*f.pvalues[k+1])
        elseif deriv ==1 # Return the derivative
            return T( dh0(t)*f.pvalues[k]/dx + dh1(t)*f.pvalues[k+1]/dx )
        else
            return zero(T) # higher order deriatives not implemented,  TODO: This should thrown an error
        end
    end
end

end # module
