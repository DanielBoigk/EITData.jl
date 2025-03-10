# EITData.jl
This is a Julia package for generating training data for the Electrical Impedance Tomography (EIT) problem. This is my Modelling seminary project.
The Data generated are meant to be used for training a neural network to solve the EIT problem. The data are generated using the finite element method (FEM) based on the Julia library Gridap.jl.

## Installation
To install the package, you can use the Julia package manager:
pkg> add https://github.com/DanielBoigk/EITData.jl
## Usage
To use the package, you can import it in your Julia script:
using EITData

## Mesh Requirements

There are some requirements regarding the mesh. 
- Please use a .msh file
- Please make sure that all boundary points 
- please make sure that the ordering of the points is such that:
    - First all the boundary points 
    - Then the other points
Reason being that I have not implemented it in such a way that it in my current implementation if I use the assemble_vector from a CellField it takes away too much performance. (Need to figure this out), thus i just want to copy the values from the dirichlet_values or the free_values and assemble the vector directly from a CellField, this however is more complicated if the values are not in the right order since since then I need to permutate the vector entries. (And I haven't implemented that yet)
- Please make sure Gmsh is installed systemwide and callable via console. 
    - This is due to a bug in gmsh.jl (maybe this is fixed right now. Need to find the Github issue again) 
