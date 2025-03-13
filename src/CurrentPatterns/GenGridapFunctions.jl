# Please generate a function that takes a vector, splits it into dirichlet_values and free_values and generates a gridap function using the constructor for a given space

function generate_gridap_function_dirichlet(vector::Array{Float64,1}, dirichlet_values::Array{Float64,1}, space::SingleFieldFESpace, N_d, m)
    gridap_function = Gridap.Fields.FEFunction(space, free_values, dirichlet_values)
    return gridap_function

end
# Please now do that for a dictionary of such vectors
function generate_gridap_functions_dirichlet(dict::Dict{Int,Array{Float64,1}}, space::SingleFieldFESpace, N_d, m)
    dict_gridap_functions = Dict{Int,Gridap.FESpace.FEFunction}()
    for (key, value) in dict
        dict_gridap_functions[key] = generate_gridap_function_dirichlet(value, dirichlet_values, space)
    end
    return dict_gridap_functions
end