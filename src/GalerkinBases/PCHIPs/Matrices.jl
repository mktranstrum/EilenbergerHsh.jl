

import SparseArrays: sparse, sparsevec

# Precompile Matrices

# Overlap matrices betewen shape functions
hh = Matrix{Rational{Int64}}(undef, 4, 4)
hdh = Matrix{Rational{Int64}}(undef, 4, 4)
dhdh = Matrix{Rational{Int64}}(undef, 4, 4)
hhh = Array{Rational{Int64}, 3}(undef, 4, 4, 4)

hh[1,:] .= [ 13//35, 11//210, 9//70, -(13//420)]
hh[2,:] .= [11//210, 1//105, 13//420, -(1//140)]
hh[3,:] .= [9//70, 13//420, 13//35, -(11//210)]
hh[4,:] .= [-(13//420), -(1//140), -(11//210), 1//105]

hdh[1,:] .= [-(1//2), 1//10, 1//2, -(1//10)]
hdh[2,:] .= [-(1//10), 0, 1//10, -(1//60)]
hdh[3,:] .= [-(1//2), -(1//10), 1//2, 1//10]
hdh[4,:] .= [1//10, 1//60, -(1//10), 0]

dhdh[1,:] .= [6/5, 1/10, -(6/5), 1/10]
dhdh[2,:] .= [1/10, 2/15, -(1/10), -(1/30)]
dhdh[3,:] .= [-(6/5), -(1/10), 6/5, -(1/10)]
dhdh[4,:] .= [1/10, -(1/30), -(1/10), 2/15]

hhh[1,1,:] .= [43//140, 97//2520, 9//140, -(43//2520)]
hhh[1,2,:] .= [97//2520, 2//315, 1//72, -(1//280)]
hhh[1,3,:] .= [9//140, 1//72, 9//140, -(1//72)]
hhh[1,4,:] .= [-(43//2520), -(1//280), -(1//72), 1//315]
hhh[2,1,:] .= [97//2520, 2//315, 1//72, -(1//280)]
hhh[2,2,:] .= [2//315, 1//840, 1//315, -(1//1260)]
hhh[2,3,:] .= [1//72, 1//315, 43//2520, -(1//280)]
hhh[2,4,:] .= [-(1//280), -(1//1260), -(1//280), 1//1260]
hhh[3,1,:] .= [9//140, 1//72, 9//140, -(1//72)]
hhh[3,2,:] .= [1//72, 1//315, 43//2520, -(1//280)]
hhh[3,3,:] .= [9//140, 43//2520, 43//140, -(97//2520)]
hhh[3,4,:] .= [-(1//72), -(1//280), -(97//2520), 2//315]
hhh[4,1,:] .= [-(43//2520), -(1//280), -(1//72), 1//315]
hhh[4,2,:] .= [-(1//280), -(1//1260), -(1//280), 1//1260]
hhh[4,3,:] .= [-(1//72), -(1//280), -(97//2520), 2//315]
hhh[4,4,:] .= [1//315, 1//1260, 2//315, -(1//840)]


function uu(knots::Knots)
    I = zeros(Int64, 12*knots.N - 8)
    J = zeros(Int64, 12*knots.N - 8)
    V = zeros(Float64, 12*knots.N - 8)

    # First knot
   
    # First row
    I[1] = 1
    J[1] = 1
    V[1] = hh[1,1]*knots.dx[1]

    I[2] = 1
    J[2] = 2
    V[2] = hh[1,2]*knots.dx[1]^2

    I[3] = 1
    J[3] = 3
    V[3] = hh[1,3]*knots.dx[1]

    I[4] = 1
    J[4] = 4
    V[4] = hh[1,4]*knots.dx[1]^2

    # Second row
    I[5] = 2
    J[5] = 1
    V[5] = hh[2,1]*knots.dx[1]^2

    I[6] = 2
    J[6] = 2
    V[6] = hh[2,2]*knots.dx[1]^3

    I[7] = 2
    J[7] = 3
    V[7] = hh[2,3]*knots.dx[1]^2

    I[8] = 2
    J[8] = 4
    V[8] = hh[2,4]*knots.dx[1]^3

    # Loop over konts 2..N-1
    for k = 2:(knots.N - 1)
        i0 = 8 + (k-2)*12

        # First row
        I[i0 + 1] = 2*k - 1
        J[i0 + 1] = 2*(k-1) - 1
        V[i0 + 1] = hh[3,1]*knots.dx[k-1]

        I[i0 + 2] = 2*k - 1
        J[i0 + 2] = 2*(k-1)
        V[i0 + 2] = hh[3,2]*knots.dx[k-1]^2
        
        I[i0 + 3] = 2*k - 1
        J[i0 + 3] = 2*k - 1
        V[i0 + 3] = hh[3,3]*knots.dx[k-1] + hh[1,1]*knots.dx[k]
        
        I[i0 + 4] = 2*k - 1
        J[i0 + 4] = 2*k
        V[i0 + 4] = hh[3,4]*knots.dx[k-1]^2 + hh[1,2]*knots.dx[k]^2

        I[i0 + 5] = 2*k - 1
        J[i0 + 5] = 2*(k+1) - 1
        V[i0 + 5] = hh[1,3]*knots.dx[k]

        I[i0 + 6] = 2*k - 1
        J[i0 + 6] = 2*(k+1)
        V[i0 + 6] = hh[1,4]*knots.dx[k]^2

        # Second row
        I[i0 + 7] = 2*k
        J[i0 + 7] = 2*(k-1) - 1
        V[i0 + 7] = hh[4,1]*knots.dx[k-1]^2

        I[i0 + 8] = 2*k
        J[i0 + 8] = 2*(k-1)
        V[i0 + 8] = hh[4,2]*knots.dx[k-1]^3
        
        I[i0 + 9] = 2*k
        J[i0 + 9] = 2*k - 1
        V[i0 + 9] = hh[4,3]*knots.dx[k-1]^2 + hh[2,1]*knots.dx[k]^2
        
        I[i0 + 10] = 2*k
        J[i0 + 10] = 2*k
        V[i0 + 10] = hh[4,4]*knots.dx[k-1]^3 + hh[2,2]*knots.dx[k]^3

        I[i0 + 11] = 2*k
        J[i0 + 11] = 2*(k+1) - 1
        V[i0 + 11] = hh[2,3]*knots.dx[k]^2

        I[i0 + 12] = 2*k
        J[i0 + 12] = 2*(k+1)
        V[i0 + 12] = hh[2,4]*knots.dx[k]^3
    end


    # Last knot
    i0 = 8 + (knots.N-2)*12

    # First Row
    I[i0 + 1] = 2*knots.N - 1
    J[i0 + 1] = 2*knots.N - 3
    V[i0 + 1] = hh[3,1]*knots.dx[knots.N-1]

    I[i0 + 2] = 2*knots.N - 1
    J[i0 + 2] = 2*knots.N - 2
    V[i0 + 2] = hh[3,2]*knots.dx[knots.N-1]^2

    I[i0 + 3] = 2*knots.N - 1
    J[i0 + 3] = 2*knots.N - 1
    V[i0 + 3] = hh[3,3]*knots.dx[knots.N-1]

    I[i0 + 4] = 2*knots.N - 1
    J[i0 + 4] = 2*knots.N 
    V[i0 + 4] = hh[3,4]*knots.dx[knots.N-1]^2

    # Second Row
    I[i0 + 5] = 2*knots.N
    J[i0 + 5] = 2*knots.N - 3
    V[i0 + 5] = hh[4,1]*knots.dx[knots.N-1]^2

    I[i0 + 6] = 2*knots.N
    J[i0 + 6] = 2*knots.N - 2
    V[i0 + 6] = hh[4,2]*knots.dx[knots.N-1]^3

    I[i0 + 7] = 2*knots.N
    J[i0 + 7] = 2*knots.N - 1
    V[i0 + 7] = hh[4,3]*knots.dx[knots.N-1]^2

    I[i0 + 8] = 2*knots.N
    J[i0 + 8] = 2*knots.N 
    V[i0 + 8] = hh[4,4]*knots.dx[knots.N-1]^3

    return sparse(I, J, V)
end

function udu(knots::Knots)
    knots.N
    I = zeros(Int64, 12*knots.N - 8)
    J = zeros(Int64, 12*knots.N - 8)
    V = zeros(Float64, 12*knots.N - 8)

    # First knot
   
    # First row
    I[1] = 1
    J[1] = 1
    V[1] = hdh[1,1]

    I[2] = 1
    J[2] = 2
    V[2] = hdh[1,2]*knots.dx[1]

    I[3] = 1
    J[3] = 3
    V[3] = hdh[1,3]

    I[4] = 1
    J[4] = 4
    V[4] = hdh[1,4]*knots.dx[1]

    # Second row
    I[5] = 2
    J[5] = 1
    V[5] = hdh[2,1]*knots.dx[1]

    I[6] = 2
    J[6] = 2
    V[6] = hdh[2,2]*knots.dx[1]^2

    I[7] = 2
    J[7] = 3
    V[7] = hdh[2,3]*knots.dx[1]

    I[8] = 2
    J[8] = 4
    V[8] = hdh[2,4]*knots.dx[1]^2

    # Loop over konts 2..N-1
    for k = 2:(knots.N - 1)
        i0 = 8 + (k-2)*12

        # First row
        I[i0 + 1] = 2*k - 1
        J[i0 + 1] = 2*(k-1) - 1
        V[i0 + 1] = hdh[3,1]

        I[i0 + 2] = 2*k - 1
        J[i0 + 2] = 2*(k-1)
        V[i0 + 2] = hdh[3,2]*knots.dx[k-1]
        
        I[i0 + 3] = 2*k - 1
        J[i0 + 3] = 2*k - 1
        V[i0 + 3] = hdh[3,3] + hdh[1,1]
        
        I[i0 + 4] = 2*k - 1
        J[i0 + 4] = 2*k
        V[i0 + 4] = hdh[3,4]*knots.dx[k-1] + hdh[1,2]*knots.dx[k]

        I[i0 + 5] = 2*k - 1
        J[i0 + 5] = 2*(k+1) - 1
        V[i0 + 5] = hdh[1,3]

        I[i0 + 6] = 2*k - 1
        J[i0 + 6] = 2*(k+1)
        V[i0 + 6] = hdh[1,4]*knots.dx[k]

        # Second row
        I[i0 + 7] = 2*k
        J[i0 + 7] = 2*(k-1) - 1
        V[i0 + 7] = hdh[4,1]*knots.dx[k-1]

        I[i0 + 8] = 2*k
        J[i0 + 8] = 2*(k-1)
        V[i0 + 8] = hdh[4,2]*knots.dx[k-1]^2
        
        I[i0 + 9] = 2*k
        J[i0 + 9] = 2*k - 1
        V[i0 + 9] = hdh[4,3]*knots.dx[k-1] + hdh[2,1]*knots.dx[k]
        
        I[i0 + 10] = 2*k
        J[i0 + 10] = 2*k
        V[i0 + 10] = hdh[4,4]*knots.dx[k-1]^2 + hdh[2,2]*knots.dx[k]^2

        I[i0 + 11] = 2*k
        J[i0 + 11] = 2*(k+1) - 1
        V[i0 + 11] = hdh[2,3]*knots.dx[k]

        I[i0 + 12] = 2*k
        J[i0 + 12] = 2*(k+1)
        V[i0 + 12] = hdh[2,4]*knots.dx[k]^2
    end


    # Last knot
    i0 = 8 + (knots.N-2)*12

    # First Row
    I[i0 + 1] = 2*knots.N - 1
    J[i0 + 1] = 2*knots.N - 3
    V[i0 + 1] = hdh[3,1]

    I[i0 + 2] = 2*knots.N - 1
    J[i0 + 2] = 2*knots.N - 2
    V[i0 + 2] = hdh[3,2]*knots.dx[knots.N-1]

    I[i0 + 3] = 2*knots.N - 1
    J[i0 + 3] = 2*knots.N - 1
    V[i0 + 3] = hdh[3,3]

    I[i0 + 4] = 2*knots.N - 1
    J[i0 + 4] = 2*knots.N 
    V[i0 + 4] = hdh[3,4]*knots.dx[knots.N-1]

    # Second Row
    I[i0 + 5] = 2*knots.N
    J[i0 + 5] = 2*knots.N - 3
    V[i0 + 5] = hdh[4,1]*knots.dx[knots.N-1]

    I[i0 + 6] = 2*knots.N
    J[i0 + 6] = 2*knots.N - 2
    V[i0 + 6] = hdh[4,2]*knots.dx[knots.N-1]^2

    I[i0 + 7] = 2*knots.N
    J[i0 + 7] = 2*knots.N - 1
    V[i0 + 7] = hdh[4,3]*knots.dx[knots.N-1]

    I[i0 + 8] = 2*knots.N
    J[i0 + 8] = 2*knots.N 
    V[i0 + 8] = hdh[4,4]*knots.dx[knots.N-1]^2

    return sparse(I, J, V)
end

function dudu(knots::Knots)
    I = zeros(Int64, 12*knots.N - 8)
    J = zeros(Int64, 12*knots.N - 8)
    V = zeros(Float64, 12*knots.N - 8)

    # First knot
   
    # First row
    I[1] = 1
    J[1] = 1
    V[1] = dhdh[1,1]/knots.dx[1]

    I[2] = 1
    J[2] = 2
    V[2] = dhdh[1,2]

    I[3] = 1
    J[3] = 3
    V[3] = dhdh[1,3]/knots.dx[1]

    I[4] = 1
    J[4] = 4
    V[4] = dhdh[1,4]

    # Second row
    I[5] = 2
    J[5] = 1
    V[5] = dhdh[2,1]

    I[6] = 2
    J[6] = 2
    V[6] = dhdh[2,2]*knots.dx[1]

    I[7] = 2
    J[7] = 3
    V[7] = dhdh[2,3]

    I[8] = 2
    J[8] = 4
    V[8] = dhdh[2,4]*knots.dx[1]

    # Loop over konts 2..N-1
    for k = 2:(knots.N - 1)
        i0 = 8 + (k-2)*12

        # First row
        I[i0 + 1] = 2*k - 1
        J[i0 + 1] = 2*(k-1) - 1
        V[i0 + 1] = dhdh[3,1]/knots.dx[k-1]

        I[i0 + 2] = 2*k - 1
        J[i0 + 2] = 2*(k-1)
        V[i0 + 2] = dhdh[3,2]
        
        I[i0 + 3] = 2*k - 1
        J[i0 + 3] = 2*k - 1
        V[i0 + 3] = dhdh[3,3]/knots.dx[k-1] + dhdh[1,1]/knots.dx[k]
        
        I[i0 + 4] = 2*k - 1
        J[i0 + 4] = 2*k
        V[i0 + 4] = dhdh[3,4] + dhdh[1,2]

        I[i0 + 5] = 2*k - 1
        J[i0 + 5] = 2*(k+1) - 1
        V[i0 + 5] = dhdh[1,3]/knots.dx[k]

        I[i0 + 6] = 2*k - 1
        J[i0 + 6] = 2*(k+1)
        V[i0 + 6] = dhdh[1,4]

        # Second row
        I[i0 + 7] = 2*k
        J[i0 + 7] = 2*(k-1) - 1
        V[i0 + 7] = dhdh[4,1]

        I[i0 + 8] = 2*k
        J[i0 + 8] = 2*(k-1)
        V[i0 + 8] = dhdh[4,2]*knots.dx[k-1]
        
        I[i0 + 9] = 2*k
        J[i0 + 9] = 2*k - 1
        V[i0 + 9] = dhdh[4,3] + dhdh[2,1]
        
        I[i0 + 10] = 2*k
        J[i0 + 10] = 2*k
        V[i0 + 10] = dhdh[4,4]*knots.dx[k-1] + dhdh[2,2]*knots.dx[k]

        I[i0 + 11] = 2*k
        J[i0 + 11] = 2*(k+1) - 1
        V[i0 + 11] = dhdh[2,3]

        I[i0 + 12] = 2*k
        J[i0 + 12] = 2*(k+1)
        V[i0 + 12] = dhdh[2,4]*knots.dx[k]
    end


    # Last knot
    i0 = 8 + (knots.N-2)*12

    # First Row
    I[i0 + 1] = 2*knots.N - 1
    J[i0 + 1] = 2*knots.N - 3
    V[i0 + 1] = dhdh[3,1]/knots.dx[knots.N-1]

    I[i0 + 2] = 2*knots.N - 1
    J[i0 + 2] = 2*knots.N - 2
    V[i0 + 2] = dhdh[3,2]

    I[i0 + 3] = 2*knots.N - 1
    J[i0 + 3] = 2*knots.N - 1
    V[i0 + 3] = dhdh[3,3]/knots.dx[knots.N-1]

    I[i0 + 4] = 2*knots.N - 1
    J[i0 + 4] = 2*knots.N 
    V[i0 + 4] = dhdh[3,4]

    # Second Row
    I[i0 + 5] = 2*knots.N
    J[i0 + 5] = 2*knots.N - 3
    V[i0 + 5] = dhdh[4,1]

    I[i0 + 6] = 2*knots.N
    J[i0 + 6] = 2*knots.N - 2
    V[i0 + 6] = dhdh[4,2]*knots.dx[knots.N-1]

    I[i0 + 7] = 2*knots.N
    J[i0 + 7] = 2*knots.N - 1
    V[i0 + 7] = dhdh[4,3]

    I[i0 + 8] = 2*knots.N
    J[i0 + 8] = 2*knots.N 
    V[i0 + 8] = dhdh[4,4]*knots.dx[knots.N-1]

    return sparse(I, J, V)
end

function uuu(knots::Knots, i)
    if i == 1
        I = zeros(Int64, 16)
        J = zeros(Int64, 16)
        V = zeros(Float64, 16)

        # First row
        I[1] = 1
        J[1] = 1
        V[1] = hhh[1,1,1]*knots.dx[1]
        
        I[2] = 1
        J[2] = 2
        V[2] = hhh[1,1,2]*knots.dx[1]^2

        I[3] = 1
        J[3] = 3
        V[3] = hhh[1,1,3]*knots.dx[1]

        I[4] = 1
        J[4] = 4
        V[4] = hhh[1,1,4]*knots.dx[1]^2

        # Second row
        I[5] = 2
        J[5] = 1
        V[5] = hhh[1,2,1]*knots.dx[1]^2
        
        I[6] = 2
        J[6] = 2
        V[6] = hhh[1,2,2]*knots.dx[1]^3

        I[7] = 2
        J[7] = 3
        V[7] = hhh[1,2,3]*knots.dx[1]^2

        I[8] = 2
        J[8] = 4
        V[8] = hhh[1,2,4]*knots.dx[1]^3

        # Third row
        I[9] = 3
        J[9] = 1
        V[9] = hhh[1,3,1]*knots.dx[1]
        
        I[10] = 3
        J[10] = 2
        V[10] = hhh[1,3,2]*knots.dx[1]^2

        I[11] = 3
        J[11] = 3
        V[11] = hhh[1,3,3]*knots.dx[1]

        I[12] = 3
        J[12] = 4
        V[12] = hhh[1,3,4]*knots.dx[1]^2

        # Fourth row
        I[13] = 4
        J[13] = 1
        V[13] = hhh[1,4,1]*knots.dx[1]^2
        
        I[14] = 4
        J[14] = 2
        V[14] = hhh[1,4,2]*knots.dx[1]^3

        I[15] = 4
        J[15] = 3
        V[15] = hhh[1,4,3]*knots.dx[1]^2

        I[16] = 4
        J[16] = 4
        V[16] = hhh[1,4,4]*knots.dx[1]^3
        return sparse(I, J, V, 2*knots.N, 2*knots.N)
    elseif i == 2
        I = zeros(Int64, 16)
        J = zeros(Int64, 16)
        V = zeros(Float64, 16)

        # First row
        I[1] = 1
        J[1] = 1
        V[1] = hhh[2,1,1]*knots.dx[1]^2
        
        I[2] = 1
        J[2] = 2
        V[2] = hhh[2,1,2]*knots.dx[1]^3

        I[3] = 1
        J[3] = 3
        V[3] = hhh[2,1,3]*knots.dx[1]^2

        I[4] = 1
        J[4] = 4
        V[4] = hhh[2,1,4]*knots.dx[1]^3

        # Second row
        I[5] = 2
        J[5] = 1
        V[5] = hhh[2,2,1]*knots.dx[1]^3
        
        I[6] = 2
        J[6] = 2
        V[6] = hhh[2,2,2]*knots.dx[1]^4

        I[7] = 2
        J[7] = 3
        V[7] = hhh[2,2,3]*knots.dx[1]^3

        I[8] = 2
        J[8] = 4
        V[8] = hhh[2,2,4]*knots.dx[1]^4

        # Third row
        I[9] = 3
        J[9] = 1
        V[9] = hhh[2,3,1]*knots.dx[1]^2
        
        I[10] = 3
        J[10] = 2
        V[10] = hhh[2,3,2]*knots.dx[1]^3

        I[11] = 3
        J[11] = 3
        V[11] = hhh[2,3,3]*knots.dx[1]^2

        I[12] = 3
        J[12] = 4
        V[12] = hhh[2,3,4]*knots.dx[1]^3

        # Fourth row
        I[13] = 4
        J[13] = 1
        V[13] = hhh[2,4,1]*knots.dx[1]^3
        
        I[14] = 4
        J[14] = 2
        V[14] = hhh[2,4,2]*knots.dx[1]^4

        I[15] = 4
        J[15] = 3
        V[15] = hhh[2,4,3]*knots.dx[1]^3

        I[16] = 4
        J[16] = 4
        V[16] = hhh[2,4,4]*knots.dx[1]^4
        return sparse(I, J, V, 2*knots.N, 2*knots.N)
    elseif i == 2*knots.N - 1
        I = zeros(Int64, 16)
        J = zeros(Int64, 16)
        V = zeros(Float64, 16)

        # First row
        I[1] = i - 2
        J[1] = i - 2
        V[1] = hhh[3,1,1]*knots.dx[end]
        
        I[2] = i - 2
        J[2] = i - 1
        V[2] = hhh[3,1,2]*knots.dx[end]^2

        I[3] = i - 2
        J[3] = i
        V[3] = hhh[3,1,3]*knots.dx[end]

        I[4] = i - 2
        J[4] = i + 1
        V[4] = hhh[3,1,4]*knots.dx[end]^2

        # Second row
        I[5] = i - 1
        J[5] = i - 2
        V[5] = hhh[3,2,1]*knots.dx[end]^2
        
        I[6] = i - 1
        J[6] = i - 1
        V[6] = hhh[3,2,2]*knots.dx[end]^3

        I[7] = i - 1
        J[7] = i
        V[7] = hhh[3,2,3]*knots.dx[end]^2

        I[8] = i - 1
        J[8] = i + 1
        V[8] = hhh[3,2,4]*knots.dx[end]^3

        # Third row
        I[9] = i
        J[9] = i - 2
        V[9] = hhh[3,3,1]*knots.dx[end]
        
        I[10] = i
        J[10] = i - 1
        V[10] = hhh[3,3,2]*knots.dx[end]^2

        I[11] = i
        J[11] = i 
        V[11] = hhh[3,3,3]*knots.dx[end]

        I[12] = i
        J[12] = i + 1
        V[12] = hhh[3,3,4]*knots.dx[end]^2

        # Fourth row
        I[13] = i + 1
        J[13] = i - 2
        V[13] = hhh[3,4,1]*knots.dx[end]^2
        
        I[14] = i + 1
        J[14] = i - 1
        V[14] = hhh[3,4,2]*knots.dx[end]^3

        I[15] = i + 1
        J[15] = i
        V[15] = hhh[3,4,3]*knots.dx[end]^2

        I[16] = i + 1
        J[16] = i + 1
        V[16] = hhh[3,4,4]*knots.dx[end]^3
        return sparse(I, J, V, 2*knots.N, 2*knots.N)
    elseif i == 2*knots.N
        I = zeros(Int64, 16)
        J = zeros(Int64, 16)
        V = zeros(Float64, 16)

        # First row
        I[1] = i - 3
        J[1] = i - 3
        V[1] = hhh[4,1,1]*knots.dx[end]^2
        
        I[2] = i - 3
        J[2] = i - 2
        V[2] = hhh[4,1,2]*knots.dx[end]^3

        I[3] = i - 3
        J[3] = i - 1
        V[3] = hhh[4,1,3]*knots.dx[end]^2

        I[4] = i - 3
        J[4] = i
        V[4] = hhh[4,1,4]*knots.dx[end]^3

        # Second row
        I[5] = i - 2
        J[5] = i - 3
        V[5] = hhh[4,2,1]*knots.dx[end]^3
        
        I[6] = i - 2
        J[6] = i - 2
        V[6] = hhh[4,2,2]*knots.dx[end]^4

        I[7] = i - 2
        J[7] = i -1
        V[7] = hhh[4,2,3]*knots.dx[end]^3

        I[8] = i - 2
        J[8] = i
        V[8] = hhh[4,2,4]*knots.dx[end]^4

        # Third row
        I[9] = i - 1
        J[9] = i - 3
        V[9] = hhh[4,3,1]*knots.dx[end]^2
        
        I[10] = i - 1
        J[10] = i - 2
        V[10] = hhh[4,3,2]*knots.dx[end]^3

        I[11] = i - 1
        J[11] = i - 1
        V[11] = hhh[4,3,3]*knots.dx[end]^2

        I[12] = i - 1
        J[12] = i
        V[12] = hhh[4,3,4]*knots.dx[end]^3

        # Fourth row
        I[13] = i
        J[13] = i - 3
        V[13] = hhh[4,4,1]*knots.dx[end]^3
        
        I[14] = i
        J[14] = i - 2
        V[14] = hhh[4,4,2]*knots.dx[end]^4

        I[15] = i
        J[15] = i - 1
        V[15] = hhh[4,4,3]*knots.dx[end]^3

        I[16] = i
        J[16] = i
        V[16] = hhh[4,4,4]*knots.dx[end]^4
        return sparse(I, J, V, 2*knots.N, 2*knots.N)
elseif i % 2 == 1
# We're on an interior knot and specifying a pvalues

    k = div(i + 1, 2) # knot
    I = zeros(Int64, 28)
    J = zeros(Int64, 28)
    V = zeros(Float64, 28)
    
    I[1] = i - 2
    J[1] = i - 2
    V[1] = hhh[3,1,1]*knots.dx[k-1]

    I[2] = i - 2
    J[2] = i - 1
    V[2] = hhh[3,1,2]*knots.dx[k-1]^2

    I[3] = i - 2
    J[3] = i
    V[3] = hhh[3,1,3]*knots.dx[k-1]

    I[4] = i - 2
    J[4] = i + 1
    V[4] = hhh[3,1,4]*knots.dx[k-1]^2

    I[5] = i - 1
    J[5] = i - 2
    V[5] = hhh[3,2,1]*knots.dx[k-1]^2

    I[6] = i - 1
    J[6] = i - 1
    V[6] = hhh[3,2,2]*knots.dx[k-1]^3

    I[7] = i - 1
    J[7] = i
    V[7] = hhh[3,2,3]*knots.dx[k-1]^2

    I[8] = i - 1
    J[8] = i + 1
    V[8] = hhh[3,2,4]*knots.dx[k-1]^3
    
    I[9] = i
    J[9] = i - 2
    V[9] = hhh[3,3,1]*knots.dx[k-1]

    I[10] = i
    J[10] = i - 1
    V[10] = hhh[3,3,2]*knots.dx[k-1]^2

    I[11] = i
    J[11] = i
    V[11] = hhh[3,3,3]*knots.dx[k-1] + hhh[1,1,1]*knots.dx[k]

    I[12] = i
    J[12] = i + 1
    V[12] = hhh[3,3,4]*knots.dx[k-1]^2 + hhh[1,1,2]*knots.dx[k]^2

    I[13] = i
    J[13] = i + 2
    V[13] = hhh[1,1,3]*knots.dx[k]

    I[14] = i
    J[14] = i + 3
    V[14] = hhh[1,1,4]*knots.dx[k]^2

    I[15] = i + 1
    J[15] = i - 2
    V[15] = hhh[3,4,1]*knots.dx[k-1]^2

    I[16] = i + 1
    J[16] = i - 1
    V[16] = hhh[3,4,2]*knots.dx[k-1]^3

    I[17] = i + 1
    J[17] = i
    V[17] = hhh[3,4,3]*knots.dx[k-1]^2 + hhh[1,2,1]*knots.dx[k]^2

    I[18] = i + 1
    J[18] = i + 1
    V[18] = hhh[3,4,4]*knots.dx[k-1]^3 + hhh[1,2,2]*knots.dx[k]^3

    I[19] = i + 1
    J[19] = i + 2
    V[19] = hhh[1,2,3]*knots.dx[k]^2

    I[20] = i + 1
    J[20] = i + 3
    V[20] = hhh[1,2,4]*knots.dx[k]^3

    I[21] = i + 2
    J[21] = i
    V[21] = hhh[1,3,1]*knots.dx[k]

    I[22] = i + 2
    J[22] = i + 1
    V[22] = hhh[1,3,2]*knots.dx[k]^2

    I[23] = i + 2
    J[23] = i + 2
    V[23] = hhh[1,3,3]*knots.dx[k]

    I[24] = i + 2
    J[24] = i + 3
    V[24] = hhh[1,3,4]*knots.dx[k]^2

    I[25] = i + 3
    J[25] = i
    V[25] = hhh[1,4,1]*knots.dx[k]^2

    I[26] = i + 3
    J[26] = i + 1
    V[26] = hhh[1,4,2]*knots.dx[k]^3

    I[27] = i + 3
    J[27] = i + 2
    V[27] = hhh[1,4,3]*knots.dx[k]^2

    I[28] = i + 3
    J[28] = i + 3
    V[28] = hhh[1,4,4]*knots.dx[k]^3
    return sparse(I, J, V, 2*knots.N, 2*knots.N)

elseif i % 2 == 0
    # We're on an interior knot and specifying an mvalues
    k = div(i, 2) # knot
    I = zeros(Int64, 28)
    J = zeros(Int64, 28)
    V = zeros(Float64, 28)
    
    I[1] = i - 3
    J[1] = i - 3
    V[1] = hhh[4,1,1]*knots.dx[k-1]^2

    I[2] = i - 3
    J[2] = i - 2
    V[2] = hhh[4,1,2]*knots.dx[k-1]^3

    I[3] = i - 3
    J[3] = i - 1
    V[3] = hhh[4,1,3]*knots.dx[k-1]^2

    I[4] = i - 3
    J[4] = i
    V[4] = hhh[4,1,4]*knots.dx[k-1]^3

    I[5] = i - 2
    J[5] = i - 3
    V[5] = hhh[4,2,1]*knots.dx[k-1]^3

    I[6] = i - 2
    J[6] = i - 2
    V[6] = hhh[4,2,2]*knots.dx[k-1]^4

    I[7] = i - 2
    J[7] = i - 1
    V[7] = hhh[4,2,3]*knots.dx[k-1]^3

    I[8] = i - 2
    J[8] = i
    V[8] = hhh[4,2,4]*knots.dx[k-1]^4
    
    I[9] = i - 1
    J[9] = i - 3
    V[9] = hhh[4,3,1]*knots.dx[k-1]^2

    I[10] = i - 1
    J[10] = i - 2
    V[10] = hhh[4,3,2]*knots.dx[k-1]^3

    I[11] = i - 1
    J[11] = i - 1
    V[11] = hhh[4,3,3]*knots.dx[k-1]^2 + hhh[2,1,1]*knots.dx[k]^2

    I[12] = i - 1
    J[12] = i
    V[12] = hhh[4,3,4]*knots.dx[k-1]^3 + hhh[2,1,2]*knots.dx[k]^3

    I[13] = i - 1
    J[13] = i + 1
    V[13] = hhh[2,1,3]*knots.dx[k]^2

    I[14] = i - 1
    J[14] = i + 2
    V[14] = hhh[2,1,4]*knots.dx[k]^3

    I[15] = i
    J[15] = i - 3
    V[15] = hhh[4,4,1]*knots.dx[k-1]^3

    I[16] = i
    J[16] = i - 2
    V[16] = hhh[4,4,2]*knots.dx[k-1]^4

    I[17] = i
    J[17] = i - 1
    V[17] = hhh[4,4,3]*knots.dx[k-1]^3 + hhh[2,2,1]*knots.dx[k]^3

    I[18] = i
    J[18] = i
    V[18] = hhh[4,4,4]*knots.dx[k-1]^4 + hhh[2,2,2]*knots.dx[k]^4

    I[19] = i
    J[19] = i + 1
    V[19] = hhh[2,2,3]*knots.dx[k]^3

    I[20] = i
    J[20] = i + 2
    V[20] = hhh[2,2,4]*knots.dx[k]^4

    I[21] = i + 1
    J[21] = i - 1
    V[21] = hhh[2,3,1]*knots.dx[k]^2

    I[22] = i + 1
    J[22] = i
    V[22] = hhh[2,3,2]*knots.dx[k]^3

    I[23] = i + 1
    J[23] = i + 1
    V[23] = hhh[2,3,3]*knots.dx[k]^2

    I[24] = i + 1
    J[24] = i + 2
    V[24] = hhh[2,3,4]*knots.dx[k]^3

    I[25] = i + 2
    J[25] = i - 1
    V[25] = hhh[2,4,1]*knots.dx[k]^3

    I[26] = i + 2
    J[26] = i
    V[26] = hhh[2,4,2]*knots.dx[k]^4

    I[27] = i + 2
    J[27] = i + 1
    V[27] = hhh[2,4,3]*knots.dx[k]^3

    I[28] = i + 2
    J[28] = i + 2
    V[28] = hhh[2,4,4]*knots.dx[k]^4
    return sparse(I, J, V, 2*knots.N, 2*knots.N)
end

end

uuu(knots::Knots) = [uuu(knots, i) for i = 1:2*knots.N]

function ufu(knots::Knots, f::AbstractVector)
    ans = f[1]*uuu(knots, 1)
    for i = 2:Base.length(f)
        ans += f[i]*uuu(knots, i)
    end
    return ans
end

δu(knots::Knots) = sparsevec([1,2*length(knots)-1], [-1.0, 1.0], 2*length(knots))

G = uu
D = udu
H = dudu
M = uuu
