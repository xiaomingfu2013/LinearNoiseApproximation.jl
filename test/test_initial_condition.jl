using LinearNoiseApproximation
using Catalyst
using Test
rn = @reaction_network model1 begin
    @parameters a b
    @species A(t) B(t)
    a, 2 * A + 3 * B --> 5 * A
    b, 3 * A + 2 * B --> 5 * B
end

u0 = [1.0, 1.0]
length_list = [5, 7, 19]
@testset for i in 0:2
    expsys = LNASystem(rn, Ω=1.0, order=i)
    @test length(expand_initial_conditions(expsys, u0)) == length_list[i+1]
end
