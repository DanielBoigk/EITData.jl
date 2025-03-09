
# Maybe try to also provide the SparseMatrixAssembler
# For whatever reasen certain functions seem to "eat" the Assembler...




struct EITOperatorData
    mesh
    γ::CellField  # Conductivity field
    Ω::Triangulation  # Domain
    dΩ::Measure  # Volume measure
    Γ::BoundaryTriangulation  # Boundary
    dΓ::Measure  # Boundary measure
    V_n::FESpace  # FE space
    U_n::FESpace  # Trial FE space
    V_d::FESpace  # FE space
    U_d::FESpace  # Trial FE space
    u_n
    v_n
    u_d
    v_d
    reffe  # Reference element
    K_d  # Dirichlet matrix
    K_n  # Neumann matrix
    # Just give conductivity as a function that is defined over [-1,1]x[-1,1]
    N_n::Int64
    N_d::Int64
    function EITOperatorData(mesh, conductivity)
        Ω = Triangulation(mesh)
        dΩ = Measure(Ω, 2)
        Γ = BoundaryTriangulation(mesh, tags= "boundary")
        dΓ = Measure(Γ, 2)
        reffe = ReferenceFE(lagrangian, Float64, 1)
        V_n = TestFESpace(mesh, reffe, conformity=:H1)
        U_n = TrialFESpace(V_n)
        V_d = TestFESpace(mesh, reffe, conformity=:H1, dirichlet_tags="boundary")
        U_d = TrialFESpace(V_d)
        γ = interpolate_everywhere(conductivity, V_n)
        a(u, v) = ∫( γₕ * ∇(v) ⋅ ∇(u) )dΩ
        assem_n = SparseMatrixAssembler(U_n, V_n)
        assem_d = SparseMatrixAssembler(U_d, V_d)
        u_n = get_trial_fe_basis(U_n)
        v_n = get_fe_basis(V_n)
        u_d = get_trial_fe_basis(U_d)
        v_d = get_fe_basis(V_n)
        matcontribs_n = a(u_n,v_n)
        matcontribs_d = a(u_d,v_d)
        matdata_n =  collect_cell_matrix(U_n,V_n,matcontribs_n)
        matdata_d =  collect_cell_matrix(U_d,V_d,matcontribs_d)
        K_n = assemble_matrix(assem_n,matdata_n)
        K_d = assemble_matrix(assem_d,matdata_d)
        N_n = num_free_dofs(V_n)
        N_d = num_free_dofs(V_d)
        new(mesh,γ,Ω,dΩ,Γ,dΓ,V_n,U_n,V_d,U_d,u_n,v_n,u_d,v_d,reffe, K_d,K_n)
    end
end