using LinearNoiseApproximation
using Test

@testset "LinearNoiseApproximation.jl" begin
    @info "Testing a simple reaction network at steady state"
    include("test_steady_state.jl")

    @info "Testing an analytic solution"
    include("test_analytic_solution.jl")

    @info "Testing the initial condition expansion"
    include("test_initial_condition.jl")
end
