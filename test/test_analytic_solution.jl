using LinearNoiseApproximation
using Test
using OrdinaryDiffEq
function sir_lna_model!(du, u, p, t)
    β, γ, N0 = p
    S, I, ϵ11, ϵ12, ϵ22 = u
    du[1] = -β * S * I / N0
    du[2] = β * S * I / N0 - γ * I
    du[3] = β * S * I / N0 - 2 * β * I * ϵ11 / N0 - 2 * β * S * ϵ12 / N0
    du[4] = -β * S * I / N0 + (-γ + β * S / N0 - β * I / N0) * ϵ12 + β * I * ϵ11 / N0 -
            β * S * ϵ22 / N0
    du[5] = γ * I + β * S * I / N0 + 2 * (β * S / N0 - γ) * ϵ22 + 2 * β * I * ϵ12 / N0
    return nothing
end
rn = @reaction_network begin
    @parameters β γ N0 # N0 is the system size
    @species S(t) I(t)
    β / N0, S + I --> 2 * I
    γ, I --> ∅
end
u0 = [999990, 10, 0, 0, 0]  # S, I, ϵ11, ϵ12, ϵ22

tspan = (0.0, 120.0)
ts = tspan[1]:1.0:tspan[2]
alg = Tsit5()
p = [
    0.4 # λ
    1 / 7 # σ
    sum(u0) # N0 system size
]
expsys = LNASystem(rn)

prob1 = ODEProblem(sir_lna_model!, u0, tspan, p)
prob2 = ODEProblem(expsys, u0[1:2], tspan, p)
sol = solve(prob1, alg; saveat=ts, abstol=1e-7, reltol=1e-7)
sol2 = solve(prob2, alg; saveat=ts, abstol=1e-7, reltol=1e-7)

@test isapprox(sol[:, end], sol2[:, end], atol=1e-5)
