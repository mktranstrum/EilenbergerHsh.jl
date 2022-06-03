abstract type DOF end

basis(F::DOF) = F.basis
dof(F::DOF) = F.dof
Base.length(F::DOF) = length(F.dof)
update!(F::DOF, values) = F.dof .= values

include(joinpath("DOFs", "PsiDOFs.jl"))
include(joinpath("DOFs", "PhiDOFs.jl"))
