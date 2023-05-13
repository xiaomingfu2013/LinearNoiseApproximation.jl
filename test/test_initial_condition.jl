using LinearNoiseApproximation
using Catalyst
using Test
const LNA = LinearNoiseApproximation

rn = @reaction_network model1 begin
    @parameters a b
    @species A(t) B(t)
    a, 2 * A + 3 * B --> 5 * A
    b, 3 * A + 2 * B --> 5 * B
end


u0 = [1.0, 1.0]
expsys = LNASystem(rn)
@test length(LNA.expand_initial_conditions(expsys, u0)) == 5
