
function default_func()
    interpolate_array_2D(gen_discrete_data_2D())
end

struct EITFullData
    #Something
    op::EITOperatorData
    data::EITData

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
            data = EITData(op, op.m-1)
        elseif modes > op.m-1
            data = EITData(op, op.m-1)
        else
            data = EITData(op,modes)
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