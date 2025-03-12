"""
    default_func()

Default conductivity function that generates discrete 2D data and interpolates it into a function.

# Returns
- A function representing the interpolated conductivity, suitable for use in `EITFullData`.
"""
function default_func()
    interpolate_array_2D(gen_discrete_data_2D())
end


"""
    struct EITFullData

A data structure that combines the operator data and boundary data for the Electrical Impedance Tomography (EIT) problem.

# Fields
- `op::EITOperatorData`: The operator data containing mesh, finite element spaces, and stiffness matrices.
- `data::EITBoundaryData`: The boundary data containing modes, singular values, basis matrices, and true solutions.

# Constructors
    EITFullData(path::String="", modes::Int64=0, cond_func::Function=default_func, args...)

Constructs an `EITFullData` instance by loading a mesh, generating or using a provided conductivity function, and creating the operator and boundary data.

# Arguments
- `path::String=""`: Path to the mesh file. If empty, a default circular mesh is created.
- `modes::Int64=0`: Number of modes to use. If `0`, uses the maximum possible modes (`m-1`).
- `cond_func::Function=default_func`: A function to generate the conductivity. Defaults to generating discrete data via `default_func`.
- `args...`: Additional arguments passed to `cond_func`.

# Alternative Constructor
    EITFullData(path::String="", modes::Int64=0, arr::Array{Float64})

Constructs an `EITFullData` instance using a provided array for conductivity, which is interpolated into a function.

# Arguments
- `path::String=""`: Path to the mesh file.
- `modes::Int64=0`: Number of modes to use.
- `arr::Array{Float64}`: A 2D or 3D array representing conductivity, which is interpolated.

# Notes
- If no mesh path is provided, a default circular mesh is generated using `create_circle_geo()` and loaded as `"circle.msh"`.
- The conductivity can be provided as a function or an array (2D or 3D), which is then interpolated using `interpolate_array_2D` or `interpolate_array_3D`.
- The number of modes is capped at `m-1`, where `m` is the number of boundary degrees of freedom from `EITOperatorData`.
- If an invalid array dimension is provided in the alternative constructor, an error message is printed, and `-1` is returned.
- This struct serves as a container for all necessary data to solve and analyze the EIT problem.
"""
struct EITFullData
    #Something
    op::EITOperatorData
    data::EITBoundaryData

    function EITFullData(path::String= "", modes::Int64=0, cond_func::Function = default_func, args...)
        #load mesh of some kind
        # maybe check if there is a path given if not
        if path==""
            #Do Something
            create_circle_geo()
            # Just load standard path.
            mesh = GmshDiscreteModel("circle.msh")
        else
            mesh = GmshDiscreteModel(path)
        end
        # if given no function look at given sigma and mode and generate interpolable sigma function

        
        conductivity = cond_func(args...)

        if typeof(conductivity) == Array
            if ndims(conductivity) == 2
                conductivity = interpolate_array_2D(conductivity)
            elseif ndims(conductivity) == 3
                conductivity = interpolate_array_3D(conductivity)
            end
        end
        
        # Now generate the Operator from that

        
        op = EITOperatorData(mesh,conductivity)


        # now generate a basis and solve
        # Do SVD and solve; Discard unnecessarry modes
        if modes <= 0 
            data = EITBoundaryData(op, op.m-1)
        elseif modes > op.m-1
            data = EITBoundaryData(op, op.m-1)
        else
            data = EITBoundaryData(op,modes)
        end
        
        new(op, data)
    end

    function EITFullData(path::String= "", modes::Int64=0, arr::Array{Float64})

        if ndims(arr) == 2
            itp = interpolate_array_2D(arr)
        elseif ndims(arr) == 3
            itp = interpolate_array_3D(arr)
        # Now generate the Operator from that
        else
            println("Wrong input dimensions, please provide a 2D or §D array")
            return -1
        end


        return EITFullData(path, modes, itp)
    end

end