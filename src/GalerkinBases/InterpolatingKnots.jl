module InterpolatingKnots
export Knots, isequel, get_interval

"""
Knots for a interpolating spline
   
"""
struct Knots
    x::AbstractVector{Float64}
    dx::AbstractVector{Float64}
    N::Integer
end

function Knots(x::AbstractVector{Float64})
    dx = diff(x)
    @assert all(dx .> 0)
    return Knots(x, dx, Base.length(x))
end

function Base.length(knots::Knots)::Int64
    return Base.length(knots.x)
end

function isequal(k1::Knots, k2::Knots)::Bool
    if Base.isequal(k1, k2)
        return true
    else
        if length(k1) == length(k2)
            return all(k1.x .== k2.x)
        else
            return false
        end
    end
end

"""
Returns the interval in knots that contains x 
"""
function get_interval(knots::Knots, x::Real)::Integer
    if x < knots.x[1] || x > knots.x[end]
        return -1
    else
        return get_interval(knots, x, 1, knots.N)
    end
end

function get_interval(knots::Knots, x::Real, kl::Integer, ku::Integer)::Integer
    if ku - kl == 1
        return kl
    else
        k = div(kl + ku, 2)
        if knots.x[k] > x
            return get_interval(knots, x, kl, k)
        else
            return get_interval(knots, x, k, ku)
        end            
    end    
end

end # module
