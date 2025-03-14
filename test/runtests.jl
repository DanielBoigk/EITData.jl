#using Pkg
#Pkg.activate(@__DIR__)   # This ensures the test environment is active

using Test
using EITData   # or "using .EITData" if it's a local module

@testset "EITData.jl" begin
    @test isdefined(EITData, :default_func)
end