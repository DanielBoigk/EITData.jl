
# Maybe try to also provide the SparseMatrixAssembler
# For whatever reason certain functions seem to "eat" the Assembler...


"""
    struct EITOperatorData

A data structure containing the operator and finite element space information for the Electrical Impedance Tomography (EIT) problem.

# Fields
- `mesh`: The computational mesh.
- `γ::CellField`: The conductivity field interpolated over the finite element space.
- `Ω::Triangulation`: The domain triangulation.
- `dΩ::Measure`: The volume measure for integration.
- `Γ::BoundaryTriangulation`: The boundary triangulation.
- `dΓ::Measure`: The boundary measure for integration.
- `V_n::FESpace`: The test finite element space for Neumann problems.
- `U_n::FESpace`: The trial finite element space for Neumann problems.
- `V_d::FESpace`: The test finite element space for Dirichlet problems.
- `U_d::FESpace`: The trial finite element space for Dirichlet problems.
- `u_n`: The trial finite element basis for Neumann problems.
- `v_n`: The test finite element basis for Neumann problems.
- `u_d`: The trial finite element basis for Dirichlet problems.
- `v_d`: The test finite element basis for Dirichlet problems.
- `reffe`: The reference finite element (Lagrangian, order 1).
- `K_d`: The stiffness matrix for Dirichlet problems.
- `K_n`: The stiffness matrix for Neumann problems.
- `N_n::Int64`: The number of free degrees of freedom for Neumann problems.
- `N_d::Int64`: The number of free degrees of freedom for Dirichlet problems.
- `m::Int64`: The number of boundary degrees of freedom (`N_n - N_d`).

# Constructor
    EITOperatorData(mesh, conductivity)

Constructs an `EITOperatorData` instance using the provided mesh and conductivity function.

# Arguments
- `mesh`: The computational mesh.
- `conductivity`: A function or array representing the conductivity, defined over the domain (e.g., [-1,1] × [-1,1]), which is interpolated onto the finite element space.

# Notes
- The constructor sets up finite element spaces, assembles stiffness matrices for both Neumann and Dirichlet problems, and computes the necessary degrees of freedom.
- The conductivity `γ` is interpolated onto the finite element space `V_n` using `interpolate_everywhere`.
- The stiffness matrices `K_n` and `K_d` are assembled using the weak form of the EIT problem: ∫(γ * ∇v ⋅ ∇u) dΩ.
- The finite element spaces use H1 conformity, with Dirichlet conditions applied on the boundary for `V_d` and `U_d`.
"""
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
    m::Int64
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
        m = N_n-N_d
        new(mesh,γ,Ω,dΩ,Γ,dΓ,V_n,U_n,V_d,U_d,u_n,v_n,u_d,v_d,reffe, K_d,K_n, N_n,N_d,m)
    end
end