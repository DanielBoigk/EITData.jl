

struct EITFullData
    #Something
    op::EITOperatorData


    function EITFullData(path::String= "", σ=1.0, mode="exp")
        #load mesh of some kind
        # maybe check if there is a path given if not
        if path==""
            #Do Something
            # Just load standard path.
            mesh = GmshDiscreteModel(path)
        else
            mesh = GmshDiscreteModel(path)
        end
        # if given no function look at given sigma and mode and generate interpolable sigma function


        conductivity = ...
        # Now generate the Operator from that

        op = EITOperatorData(mesh,conductivity)


        # now generate a basis and solve
        # Do SVD and solve; Discard unnecessarry modes
        new(...)
    end

end