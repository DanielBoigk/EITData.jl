
"""
    make_basis_matrices(K_n, N_n, N_d, useSVD::Bool=true)

Generate basis matrices for the Electrical Impedance Tomography (EIT) problem by creating orthonormal, mean-zero basis vectors and applying necessary transformations.

# Arguments
- `K_n`: The stiffness matrix (or system matrix) used in the EIT problem.
- `N_n::Int`: The total number of nodes or elements in the system.
- `N_d::Int`: The number of interior degrees of freedom.
- `useSVD::Bool=true`: If true, applies Singular Value Decomposition (SVD) to orthogonalize the basis matrices; otherwise, returns preliminary matrices without orthogonalization.

# Returns
- If `useSVD` is `true`:
  - `G_prelim::Matrix{Float64}`: The orthogonalized basis matrix for the boundary.
  - `F_prelim::Matrix{Float64}`: The orthogonalized basis matrix for the interior.
  - `S::Vector{Float64}`: The singular values from the SVD.
- If `useSVD` is `false`:
  - `G_prelim::Matrix{Float64}`: The preliminary basis matrix for the boundary.
  - `F_prelim::Matrix{Float64}`: The preliminary basis matrix for the interior.
  - `nothing`: No singular values are computed.

# Notes
- The parameter `m = N_n - N_d` represents the number of boundary degrees of freedom.
- Preliminary basis vectors are generated using `generate_basis_vectors(N_n, m-1)`.
- The system `K_n * U_n = G` is solved using LU factorization of `K_n`.
- The function ensures zero-mean vectors over boundary points using `subtract_column_mean!`.
- When `useSVD` is `true`, SVD is applied to `F_prelim * G_prelim'`, followed by QR decomposition to ensure orthonormality of the basis matrices.
- The function assumes mesh points are ordered with the first `m` points corresponding to the boundary; incorrect ordering may lead to erroneous results.
- This function is critical for creating basis matrices representing boundary and interior potentials in EIT.
"""

function make_basis_matrices(K,N_n,N_d, useSVD::Bool = true)
    m = N_n-N_d
    G = generate_basis_vectors(N_n, m-1)
    K_factorized = lu(K)
    U_n = K_factorized \ G

    U_n = subtract_column_mean!(U_n,m)
    # This will lead to incorrect matrix if the mesh points are not ordered as described.
    G_prelim = zeros(m,m)
    F_prelim = zeros(m,m)
    G_prelim[:,m] = (1/sqrt(m))*ones(m)
    G_prelim[:,1:m-1] = G[1:m,1:m-1]
    F_prelim[:,1:m-1] = U_n[1:m,1:m-1]
    if useSVD
        U, S, V = svd(F_prelim * G_prelim')
        G_prelim[:,1:m-1] =  V[:,1:m-1]
        F_prelim[:,1:m-1] = U[:,1:m-1]
        F_prelim[:, 1:m-1] = qr(F_prelim).Q[:, 1:m-1]
        G_prelim[:, 1:m-1] = qr(G_prelim).Q[:, 1:m-1]
        F_prelim[:, 1:m-1] = subtract_column_mean!(F_prelim[:, 1:m-1],m)
        G_prelim[:, 1:m-1] = subtract_column_mean!(G_prelim[:, 1:m-1],m)
        #please multiply columns of F_prelim by the corresponding singular values from S
        F_prelim[:, 1:m-1] = F_prelim[:, 1:m-1] * Diagonal(S[1:m-1]) 
        return G_prelim, F_prelim, S
    else
        # it should not enter this branch if using the standard constructor for EITBoundaryData, since it useSVD is true by default
        return G_prelim, F_prelim, nothing
    end
end


"""
    make_force_vectors(M_short::Array{Float64,2}, N::Int)

Create a dictionary of force vectors by extending a short matrix `M_short` to full size and storing each column as a separate vector.

# Arguments
- `M_short::Array{Float64,2}`: A matrix of size `(m, modes)`, where `m` is the number of boundary points and `modes` is the number of force vectors.
- `N::Int`: The total number of nodes or elements in the system.

# Returns
- `M_dict::Dict{Int, Array{Float64,2}}`: A dictionary where each key `i` (from 1 to `modes`) maps to a force vector of size `(N, 1)`.

# Notes
- The number of modes is determined by `modes = size(M_short, 2)`.
- A full-sized matrix `M` of size `(N, modes)` is created, with `M_short` placed in the first `m` rows and the remaining rows filled with zeros.
- Each column of `M` is extracted and stored in a dictionary with integer keys from 1 to `modes`.
- This function is useful in EIT simulations for organizing force vectors, each representing a distinct boundary condition or mode.
"""
function make_force_vectors(M_short::Array{Float64,2}, N::Int)
    m = size(M_short,2)
    M = zeros(Float64, N, m)
    M[1:m, :] = M_short
    # please also put every vector in a dictionary with labels 1, ...,m
    M_dict = Dict{Int,Array{Float64,2}}()
    for i in 1:m
        M_dict[i] = M[:, i:i]
    end
    M_dict
end

"""
    calc_true_solutions(K::Array{Float64,2}, G::Array{Float64,2}, V::Gridap.FESpaces.FESpace, m, modes)

Calculate the true solutions for the EIT problem by solving the linear system for each mode and storing the results as finite element functions.

# Arguments
- `K::Array{Float64,2}`: The stiffness matrix (or system matrix) for the EIT problem.
- `G::Array{Float64,2}`: The basis matrix for the boundary conditions.
- `V::Gridap.FESpaces.FESpace`: The finite element space used for the solution.
- `m::Int`: The number of boundary points.
- `modes::Int`: The number of modes (basis vectors) to compute solutions for.

# Returns
- `U_dict::Dict{Int, Gridap.FESpaces.FEFunction}`: A dictionary where each key `i` (from 1 to `modes`) maps to a finite element function representing the solution for that mode.

# Notes
- The system `K * U_n = G_full` is solved using LU factorization of `K` for each mode.
- `G_full` is a matrix of size `(size(K, 1), modes)`, with the first `m` rows populated by `G[:, 1:modes]` and the rest filled with zeros.
- Solutions `U_n` are adjusted to have zero mean over the boundary points using `subtract_column_mean!`.
- Each solution vector is converted to a `Gridap.FESpaces.FEFunction` for finite element computations.
- This function is key to generating training data in EIT by solving the forward problem for various boundary conditions.
"""
function calc_true_solutions(K::Array{Float64,2}, G::Array{Float64,2},V::Gridap.FESpaces.FESpace, m, modes)
    K_LU = lu(K)
    G_full = zeros(Float64, size(K,1), modes)
    G_full[1:m, :] = G[:,1:modes]
    U_n = K_LU \ G
    U_n = subtract_column_mean!(U,m)
    U_dict = Dict{Int,Gridap.FESpaces.FEFunction}()
    for i in 1:modes
        U_dict[i] = Gridap.FESpaces.FEFunction(V, U_n[:,i])
    end
    U_dict
end