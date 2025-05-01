
# This contains all the information a potential neural net needs to know for training
struct MinimalPair
    dim::Int
    # coordinates of the mesh
    boundary_points::Array{Float64,2}
    interior_points::Array{Float64,2}
    γ_boundary::Array{Float64,1}
    γ_interior::Array{Float64,1}
    g_boundary::Array{Float64,1}
    f_boundary::Array{Float64,1}
    u_f::Array{Float64,1}
    u_g::Array{Float64,1}
    gradient_boundary::Array{Float64,2}
    gradient_interior::Array{Float64,2}
    error_L2::Float64
    error_Linf::Float64
    error_wasserstein::Float64

    #function MinimalPair(data::EITFullData)

    #end
end
