module EITData

    using Random, Statistics, LinearAlgebra, Images, Interpolations
    using Gridap, GridapGmsh, SparseArrays
    using Gridap.FESpaces
    using Gridap.Geometry
    using Gridap.ReferenceFEs
    using LineSearches: BackTracking

    #using PyCall
    #const mymodule = pyimport("../python/CalderonFEM.py")
    using SciPy
    # Non necessity to do this:
    #const MeshIO = pyimport("meshio")    
    #const np = pyimport("numpy")


    # This generates the conductivity data using gaussian filters and returns either an array or function that uses linear interpolation.
    include("GenerateConductivity/GenerateConductivity.jl")

    
    # include("Interpolation/Interpolation.jl")
    
    #include("CalculateGradient/CalculateGradient.jl")

    #include("GenerateFullData/GenerateFullData.jl")

    #include("Optimizer/Optimizer.jl")

    include("GenerateMesh/CreateGeo.jl")

    include("GenerateFullData/GenerateFullData.jl")



end # module EITData
