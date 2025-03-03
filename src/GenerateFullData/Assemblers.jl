export assemble_matrices
export construct_DtN_vector
export construct_NtD_vector

function assemble_matrices(γ, U_n, V_n, v_n, u_n, U_d, V_d, v_d, u_d, assem_n, assem_d)
    a(u, v) = ∫( γ * ∇(v) ⋅ ∇(u) )dΩ
    matcontribs_n = a(u_n, v_n)
    matdata_n = collect_cell_matrix(U_n, V_n, matcontribs_n)
    K_n = assemble_matrix(assem_n, matdata_n)
    matcontribs_d = a(u_d, v_d)
    matdata_d = collect_cell_matrix(U_d, V_d, matcontribs_d)
    K_d = assemble_matrix(assem_d, matdata_d)
    return K_n, K_d
end

function assemble_matrices(γ, U_n, V_n, v_n, u_n, U_d, V_d, v_d, u_d)
    a(u, v) = ∫( γ * ∇(v) ⋅ ∇(u) )dΩ
    assem_n = SparseMatrixAssembler(U_n, V_n)
    matcontribs_n = a(u_n, v_n)
    matdata_n = collect_cell_matrix(U_n, V_n, matcontribs_n)
    K_n = assemble_matrix(assem_n, matdata_n)
    assem_d = SparseMatrixAssembler(U_d, V_d)
    matcontribs_d = a(u_d, v_d)
    matdata_d = collect_cell_matrix(U_d, V_d, matcontribs_d)
    K_d = assemble_matrix(assem_d, matdata_d)
    return K_n, K_d
end


function construct_NtD_vector(V,v_n,g,assem,dΓ)
    l(v) = ∫( g * v )dΓ  # Neumann condition
    veccontribs = l(v_n)
    vecdata =  collect_cell_vector(V,veccontribs)
    vec = assemble_vector(assem,vecdata)
    # Ensure mean zero
    vec_mean = mean(vec)
    vec .-=vec_mean
    vec
end

function construct_DtN_vector(V, veccontribs, f::FEFunction)
    Ui = TrialFESpace(V, f)
    assem = SparseMatrixAssembler(Ui, Ui)
    vecdata =  collect_cell_vector(Ui,veccontribs)
    assemble_vector(assem,vecdata)
end

function construct_DtN_vector(V, veccontribs, u_f::AbstractArray)
    f = FEFunction(V,u_f)
    Ui = TrialFESpace(V, f)
    assem = SparseMatrixAssembler(Ui, Ui)
    vecdata =  collect_cell_vector(Ui,veccontribs)
    assemble_vector(assem,vecdata)
end

function construct_DtN_vector(V, f)
    U = TrialFESpace(V, f)
    a(u, v) = ∫( γ * ∇(v) ⋅ ∇(u) )dΩ
    l(v) = 0
    AffineFEOperator(a,l,U, V).op.vector
end

