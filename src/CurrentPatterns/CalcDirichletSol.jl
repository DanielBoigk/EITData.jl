
function subtract_column_mean!(A::Matrix, n::Int)
    # Calculate the mean of the first n entries of each column
    col_means = mean(A[1:n, :], dims=1)

    # Subtract the mean from each element of the column
    A .-= col_means
    
    return A
end

function make_basis_matrices(K_n,N_n,N_d, useSVD::bool = true)
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
        return G_prelim, F_prelim, S
    else
        # it should not enter this branch if using the standard constructor for EITBoundaryData, since it useSVD is true by default
        return G_prelim, F_prelim, nothing
    end
end

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
