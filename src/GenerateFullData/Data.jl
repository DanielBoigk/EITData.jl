export FullGridapData


# generates some preliminary orthonormal basis
function mean_zero_basis(n::Int, m::Int)
    out = zeros(n)  
    norm = inv(sqrt((m + 1) * m))
    @inbounds for i in 1:m
        out[i] = norm 
    end
    out[m + 1] = -m * norm 
    out
end
# 
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

struct FullGridapData
    mesh            # Specify the type here if you know it, e.g., Mesh
    reffe           # Element type
    Ω
    dΩ
    Γ
    dΓ

    boundary_tags

    f_dim::Int
    g_dim::Int

    U_n
    V_n
    U_d
    V_d
    K_n
    K_d

    u_n
    v_n
    u_d
    v_d

    γ               # interpolable version of γ
    γ_vec::Vector{Any}

    no_pairs::Int

    f_vectors::Vector{Any}
    g_vectors::Vector{Any}  # This should be an orthonormal basis and ordered after singular values

    u_funcs::Vector{Any}
    u_vecs::Vector{Any}

    singular_values::Vector{Any}
    
end