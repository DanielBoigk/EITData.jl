
"""
    struct EITBoundaryData

A data structure containing the boundary data for the Electrical Impedance Tomography (EIT) problem, including orthonormal basis matrices, force vectors, and true solutions.

# Fields
- `modes::Int64`: The number of modes used to generate the boundary data.
- `singular_values::Array{Float64,1}`: The singular values corresponding to the modes.
- `G::Array{Float64,2}`: Orthonormal matrix containing Neumann boundary data for the EIT problem.
- `F::Array{Float64,2}`: Matrix containing Dirichlet boundary data for the EIT problem.
- `g_vecs::Dict{Int,Array{Float64,2}}`: Dictionary of force vectors for Neumann boundary data, keyed by mode index.
- `f_vecs::Dict{Int,Array{Float64,2}}`: Dictionary of force vectors for Dirichlet boundary data, keyed by mode index.
- `u::Dict{Int,Gridap.FESpaces.SingleFieldFEFunction}`: Dictionary of true solutions (finite element functions) for each mode.

# Constructor
    EITBoundaryData(op, modes)

Constructs an `EITBoundaryData` instance using the provided operator data and number of modes.

# Arguments
- `op`: An instance of `EITOperatorData` containing the necessary operator and mesh information.
- `modes::Int64`: The number of modes to use for generating the boundary data.

# Notes
- The constructor computes orthonormal basis matrices `G` and `F` using `make_basis_matrices`, generates force vectors using `make_force_vectors`, and calculates true solutions using `calc_true_solutions`.
- The matrices `G` and `F` are orthonormal and represent the Neumann and Dirichlet boundary conditions, respectively.
- The true solutions `u` are stored as finite element functions for each mode.
- Commented-out fields `g_func` and `f_func` suggest potential future extensions for boundary data as functions, but they are not currently implemented.
"""
struct EITBoundaryData
    
    
    # This is the number of modes used to generate the boundary data.
    modes::Int64
    # The singular values for the modes
    singular_values::Array{Float64,1}
    # This matrix contains the neumann boundary data for the EIT problem.
    # This matrix is orthonormal.
    G::Array{Float64,2}
    # Values of the dirichlet boundary data for the EIT problem.
    F::Array{Float64,2}

    # force vectors for the neumann boundary data.
    g_vecs::Dict{Int,Array{Float64,2}}
    # force vectors for the dirichlet boundary data.
    f_vecs::Dict{Int,Array{Float64,2}}
    
    
    # the true solutions:
    u::Dict{Int,Gridap.FESpaces.SingleFieldFEFunction}

    # the values on the boundary as a function if that's possible:
    #g_func::Dict{Int,Gridap.FESpaces.SingleFieldFEFunction}
    #f_func::Dict{Int,Gridap.FESpaces.SingleFieldFEFunction}


    function EITBoundaryData(op, modes)
        
        N_n = op.N_n
        N_d = op.N_d
        G, F, singular_values = make_basis_matrices(op.K_n,N_n,N_d)

        g_vecs = make_force_vectors(G[:,1:modes],N_n)
        f_vecs = make_force_vectors(F[:,1:modes],N_d)

        u = calc_true_solutions(op.K_d, F[:,1:modes],op.V_d,op.m_d,modes)


        new(modes,singular_values,G,F,g_vecs,f_vecs,u 
        #,g_func,f_func
        )
    end
end