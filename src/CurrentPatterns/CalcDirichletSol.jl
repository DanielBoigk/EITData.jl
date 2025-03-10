
function subtract_column_mean!(A::Matrix, n::Int)
    # Calculate the mean of the first n entries of each column
    col_means = mean(A[1:n, :], dims=1)

    # Subtract the mean from each element of the column
    A .-= col_means
    
    return A
end

function make_basis_matrices(V_n,V_d,K_n,N_n,N_d, useSVD::bool = true)
    m = N_n-N_d
    G = generate_basis_vectors(N_n, m-1)
    K_factorized = lu(K)
    U_n = K_factorized \ G

    U_n = subtract_column_mean!(U_f,m)

    if useSVD
        #Do somehing
    end

    
end


