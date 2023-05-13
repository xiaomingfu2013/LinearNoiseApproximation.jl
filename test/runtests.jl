using LinearNoiseApproximation
using Test

@testset "LinearNoiseApproximation.jl" begin
    @info "Testing a simple reaction network"
    include("test_example.jl")

    @info "Testing the initial condition expansion"
    include("test_initial_condition.jl")
end
