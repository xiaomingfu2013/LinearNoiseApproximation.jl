using Test
using LinearNoiseApproximation
using OrdinaryDiffEq
using LinearAlgebra
const LNA = LinearNoiseApproximation

rn = @reaction_network TwoStage begin
    @parameters v0 v1 d0 d1 Ω
    @species M(t) P(t)
    v0 * Ω, ∅ --> M
    v1, M --> M + P
    d0, M --> ∅
    d1, P --> ∅
end
@variables v0, v1, d0, d1, Ω
lnasys = LNASystem(rn)
rates = [v0 => 4.0, v1 => 10.0, d0 => 1.0, d1 => 1.0, Ω => 2.0]
tspan = (0.0, 100.0)
u0 = [1.0, 1.0]
prob = ODEProblem(lnasys, u0, tspan, rates)
sol = solve(prob, Vern7(); abstol=1e-7, save_everystep=false)

soltf = sol(tspan[2])
equilibre = soltf[1:2]

A, BB = LNA.get_val_matrix_at_equilirium(rn, rates, equilibre)
C = lyap(A, BB)
C_ = soltf[3:end]
@test vcat(C...)[[1, 3, 4]] ≈ C_ atol = 1e-2
