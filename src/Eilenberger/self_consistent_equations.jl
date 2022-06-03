function self_consistent_residuals!(r::Vector, Ψ::ΨDOF{PCHIP}, parameters::EilenbergerParameters)
    r .= 0
    Φs = [ΦDOF(basis(Ψ)) for _ = 1:Threads.nthreads()]
    MΔ_mat = MΔ(Ψ)
    MAy_mat = MAy(Ψ)

    for job in [ Threads.@spawn(fermi_surface_residuals!(n,i,j,parameters,MΔ_mat, MAy_mat, Ψ, Φs)) for (n,i,j) = parameters.nij]
        wait(job)
        r .+= fetch(job)
    end
    r .*= 2π*parameters.T
    b = basis(Ψ)
    r .+= [b.G*Δhat(Ψ)*log(parameters.t);
          (parameters.κ₀^2/3) * (b.H[Ψ.Aydof, :]*Ayhat(Ψ)  + b.u0[Ψ.Aydof]*parameters.Ha)]
    nothing        
end

function fermi_surface_residuals!(n,i,j, parameters, MΔ_mat, MAy_mat, Ψ, Φs)
    Φ = Φs[Threads.threadid()]
    b = basis(Φ)
    ω, ny, nz, ρ0, ρ1, ρ2, ρ3 = ωnρ(parameters, n, i, j)
    solve_eilenberger!(Φ, MΔ_mat, MAy_mat, ω, ny, nz, parameters.ΔT)        
    return [b.G*Δhat(Ψ)*(ρ0 + ρ1 + ρ2 + ρ3)/ω - b.G*( (ρ0 + ρ3)*real.(fhat(Φ)) + (ρ1 + ρ2)*real.(fbhat(Φ)));
            -ny*b.G[Ψ.Aydof, :]*( (ρ0 + ρ1 + ρ2 + ρ3)*imag.(ghat(Φ)))]
end    

function residuals!(r, dofs, Ψ, parameters)
    update!(Ψ, dofs)
    self_consistent_residuals!(r, Ψ, parameters)
    nothing
end

function residuals(dofs, Ψ, parameters)
    r = zero(dofs)
    residuals!(r, dofs, Ψ, parameters)
    return r
end

function solve!(Ψ, parameters)
    @time result = nlsolve((r, dofs) -> residuals!(r, dofs, Ψ, parameters), dof(Ψ))
    update!(Ψ, result.zero)
    nothing
end
