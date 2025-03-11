module EITData

    using Random, Statistics, LinearAlgebra, Images, Interpolations
    using Gridap, GridapGmsh, SparseArrays
    using Gridap.FESpaces
    using Gridap.Geometry
    using Gridap.ReferenceFEs
    using LineSearches: BackTracking
    using Serialization
    #using PyCall
    #const mymodule = pyimport("../python/CalderonFEM.py")
    using SciPy
    # Non necessity to do this:
    #const MeshIO = pyimport("meshio")    
    #const np = pyimport("numpy")


    # This generates the conductivity data using gaussian filters and returns either an array or function that uses linear interpolation.
    include("GenerateConductivity/GenerateConductivity.jl")

    # This generates the mesh and the boundary conditions.
    include("CreateMesh/CreateMesh.jl")

    # This generates the full data for the EIT problem.
    include("GenerateFullData/GenerateFullData.jl")

    # This generates the basis vectors for the current patterns.
    include("CurrentPatterns/CurrentPatterns.jl")

    include("Export/Export.jl")

end # module EITData
