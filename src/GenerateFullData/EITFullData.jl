"""
    default_func()

Default conductivity function that generates discrete 2D data and interpolates it into a function.

# Returns
- A function representing the interpolated conductivity, suitable for use in `EITFullData`.
"""
function default_func()
    interpolate_array_2D(gen_discrete_data_2D())
end


# Just put in a function that can be evaluated at any point in the domain[-1,1]x[-1,1]. Ignore the rest.
struct EITFullData
    op::EITOperatorData
    data::EITBoundaryData

    function EITFullData(operator, data)
        return new(operator, data)
        
    end

    function EITFullData(cond_func, mesh, modes::Int64=0)
        op = EITOperatorData(mesh, cond_func)
        if modes <= 0 
            data = EITBoundaryData(op, op.m-1)
        elseif modes > op.m-1
            data = EITBoundaryData(op, op.m-1)
        else
            data = EITBoundaryData(op, modes)
        end
        new(op, data)
    end

    function EITFullData(cond_func = default_func(), path::String= "", modes::Int64=0)
        # load mesh of some kind
        # maybe check if there is a path given if not
        if path == ""
            # Do Something
            #Please test if circle.msh exists:
            if !isfile("circle.msh")
                create_circle_geo()
            end
            # Just load standard path.
            mesh = GmshDiscreteModel("circle.msh")
        else
            mesh = GmshDiscreteModel(path)
        end
        #=
        if typeof(cond_func) == Array
            if ndims(cond_func) == 2
                cond_func = interpolate_array_2D(cond_func)
            elseif ndims(cond_func) == 3
                cond_func = interpolate_array_3D(cond_func)
            else
                println("Conductivity function is not a function")
            return -1
           end
        end

        # Call the conductivity function without additional arguments
        conductivity = cond_func()
        
        if typeof(conductivity) == Array
            if ndims(conductivity) == 2
                conductivity = interpolate_array_2D(conductivity)
            elseif ndims(conductivity) == 3
                conductivity = interpolate_array_3D(conductivity)
            end
        end
        =#
        # Now generate the Operator from that
        op = EITOperatorData(mesh, cond_func)
        
        # now generate a basis and solve
        # Do SVD and solve; Discard unnecessarry modes
        if modes <= 0 
            data = EITBoundaryData(op, op.m-1)
        elseif modes > op.m-1
            data = EITBoundaryData(op, op.m-1)
        else
            data = EITBoundaryData(op, modes)
        end
        
        new(op, data)
    end


end

using Interpolations

function FromImageData(img)
    img_float = Float64.(img)
    # Ensure positive values by adding 1e-6 to all zero values
    img_float[img_float .== 0] .+= 1e-6
    # Interpolate the image
    img_interp = Interpolations.interpolate(img_float, BSpline(Linear()))
    # get the size of the image in each dimension
    size_x = size(img, 1)
    size_y = size(img, 2)
    domain = (-1,1,-1,1)
    partition = (size_x,size_y)
    mesh = CartesianDiscreteModel(domain,partition)
    EITFullData(img_interp,mesh)
end
