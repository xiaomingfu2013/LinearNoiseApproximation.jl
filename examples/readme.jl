# # Tutorial: The complete example from the Readme
#
# This example shows the complete code in order to run the code snippets from
# the Readme.

# The imports where omitted for brevity in the Readme
using LinearNoiseApproximation
using DifferentialEquations
using Catalyst
using Plots

# Define the reaction network with the Catalyst.jl package
rn = @reaction_network TwoStage begin
    @parameters  v0 v1 d0 d1 Ω
    @species M(t) P(t)
    v0*Ω, ∅ --> M
    v1, M --> M + P
    d0, M --> ∅
    d1, P --> ∅
end

# Automatically create the linear noise approximation
LNAsys = LNASystem(rn)

# Solve the LNA system with the DifferentialEquations.jl package
@variables v0, v1, d0, d1, Ω
rates = [v0=>4.0,v1=>10.0,d0=>1.0,d1=>1.0,Ω=>5.0]
tspan = (0.0, 20.0)
u0 =[0, 0] # initial condition for M and P
prob = ODEProblem(LNAsys, u0, tspan, rates)
sol = solve(prob,Vern7(),abstol=1e-7, saveat =1.0)

# Find the indices of the means and variances of the system
mean_idxs, var_idxs = find_states_cov_number([1, 2], LNAsys)

# Do some simple plotting of the numerical solutions, with the solid line
# showing the mean values and the shaded areas showing the standard deviation
# around the mean.
plt = Plots.plot(; size=(800, 350))
for (i, (mean_idx, var_idx)) in enumerate(zip(mean_idxs, var_idxs))
    mean = sol[mean_idx, :]
    var = sol[var_idx, :]

    Plots.plot!(
        sol.t,
        mean,
        ribbon=sqrt.(var),
    )
end
plt
