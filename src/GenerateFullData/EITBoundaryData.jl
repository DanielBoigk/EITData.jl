
# This is the data structure for the EIT boundary data.
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
    g_func::Dict{Int,Gridap.FESpaces.SingleFieldFEFunction}
    f_func::Dict{Int,Gridap.FESpaces.SingleFieldFEFunction}


    function EITBoundaryData(op, modes)
        
        N_n = op.N_n
        N_d = op.N_d
        G, F, singular_values = make_basis_matrices(op.K_n,N_n,N_d)

        g_vecs = make_force_vectors(G[:,1:modes],N_n)
        f_vecs = make_force_vectors(F[:,1:modes],N_d)
        m = op.m-1
        G_prelim = generate_basis_vectors(N_n,m)
        
        K_factorized = lu(op.K_n)

        U_f_prelim = K_factorized \ G_prelim
        U_f_prelim = subtract_column_mean!(U_f_prelim,m)
        G_small = G_prelim[:,1:m]
        U_f_small = U_f_prelim[1:m,:]




        new(modes,singular_values,G,F,g_vecs,f_vecs,u,g_func,f_func)
    end
end