# EITData.jl
This is a Julia package for generating training data for the Electrical Impedance Tomography (EIT) problem. This is my Master Thesis project.
The Data generated are meant to be used for training a neural network to solve the EIT problem. The data are generated using the finite element method (FEM) based on the Julia library Gridap.jl.

## Installation
To install the package, you can use the Julia package manager:
pkg> add https://github.com/DanielBoigk/EITData.jl
## Usage
To use the package, you can import it in your Julia script:
using EITData

