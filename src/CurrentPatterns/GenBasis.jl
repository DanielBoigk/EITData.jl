export mean_zero_basis, generate_basis_vectors


# Generates Basis vectors that are orthonormal and mean 0
function mean_zero_basis(n::Int, m::Int)
    out = zeros(n)  
    norm = inv(sqrt((m + 1) * m))
    @inbounds for i in 1:m
        out[i] = norm 
    end
    out[m + 1] = -m * norm 
    out
end

function generate_basis_vectors(n::Int, i::Int)
    basis_matrix = zeros(n, i)
    for m in 1:i
        basis_matrix[:, m] = mean_zero_basis(n, m)
    end
    basis_matrix
end

function subtract_column_mean!(A::Matrix, n::Int)
    # Calculate the mean of the first n entries of each column
    col_means = mean(A[1:n, :], dims=1)

    # Subtract the mean from each element of the column
    A .-= col_means
    
    return A
end