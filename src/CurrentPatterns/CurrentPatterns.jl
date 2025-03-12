include("GenBasis.jl")
export mean_zero_basis, generate_basis_vectors

include("CalcDirichletSol.jl")

export make_basis_matrices, make_force_vectors, calc_true_solutions
include("GenGridapFunctions.jl")
export generate_gridap_function_dirichlet, generate_gridap_functions_dirichlet