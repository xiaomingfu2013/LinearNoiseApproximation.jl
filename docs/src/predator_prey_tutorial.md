```@meta
EditURL = "../../examples/predator_prey_tutorial.jl"
```

# Tutorial: Using LinearNoiseApproximation.jl to solve a Lotka-Volterra model

This example demonstrates how to use `LinearNoiseApproximation.jl` package to solve a Lotka-Volterra model using the Linear Noise Approximation (LNA) and compare it to the stochastic
trajectories obtained using the Gillespie algorithm.

````@example predator_prey_tutorial
using LinearNoiseApproximation
using DifferentialEquations
using Catalyst
using JumpProcesses
using Plots
````

## Define the Lotka-Volterra model
The Lotka-Volterra model describes the dynamics of predator-prey interactions in an ecosystem. It assumes that the prey population $U$ grows at a rate represented by the parameter $\alpha$ in the absence of predators, but decreases as they are consumed by the predator population $V$ at a rate determined by the interaction strength parameter, $\beta$. The predator population decreases in size if they cannot find enough prey to consume, which is represented by the mortality rate parameter, $\delta$.
The corresponding ODE system reads
```math
\begin{equation}
	\begin{aligned}
		\frac{\mathrm{d}U}{\mathrm{d} t} & = \alpha U - \beta U V,\\
		\frac{\mathrm{d}V}{\mathrm{d} t} & = \beta U V  - \delta V,
	\end{aligned}\quad t\in (0, t_{\text{end}}).
\end{equation}
```

````@example predator_prey_tutorial
rn = @reaction_network begin
    @parameters α β δ
    @species U(t) V(t)
    α, U --> 2*U
    β, U + V --> 2*V
    δ, V --> 0
end
````

## Convert the model to a JumpSystem for the Gillespie algorithm
Using `Catalyst.jl` we can convert `ReactionSystem` to a `JumpSystem` which can be used to simulate the stochastic trajectories using the Gillespie algorithm.

````@example predator_prey_tutorial
jumpsys = convert(JumpSystem, rn)

#define the initial conditions, parameters, and time span
u0 = [120.0, 140.0]
ps = [0.8, 0.005, 0.4]
duration = 30.0
tspan = (0.0, duration)

#solve the model using the Gillespie algorithm
dprob = DiscreteProblem(jumpsys, u0, tspan, ps)
jprob = JumpProblem(jumpsys, dprob, Direct(), save_positions=(false, false))
jsol = solve(jprob, SSAStepper(), saveat=0.5, seed=2)

# solve the model using the Linear Noise Approximation (LNA)
lna = LNASystem(rn)
alg = Vern7()
prob = ODEProblem(lna, u0, tspan, ps)
sol = solve(prob, alg, saveat=0.5, seed=2)

#find the indices of the mean and variance states
mean_idxs, var_idxs = find_states_cov_number([1, 2], lna)


#plot the results
plt = Plots.plot(; size=(800, 350))
color_palette=["#1f78b4" "#ff7f00"]
label = ["U" "V"]

for (i, (mean_idx, var_idx)) in enumerate(zip(mean_idxs, var_idxs))
    mean = sol[mean_idx, :]
    var = sol[var_idx, :]

    Plots.plot!(
        jsol.t,
        jsol[mean_idx, :];
        label=string("SSA: ", label[i]),
        xlabel="Time",
        ylabel="Numbers",
        color=color_palette[i],
        legend=:topleft,
        st=:scatter,
        alpha=0.8,
        left_margin=5Plots.mm,
        bottom_margin=5Plots.mm,
        markerstrokewidth=0,
        xlim=(0, duration + 1),
        xticks=0:5:duration + 1,
    )
    Plots.plot!(
        sol.t,
        mean;
        label=string("LNA: ", label[i]),
        ribbon=sqrt.(var),
        ribbonalpha=0.3,
        color=color_palette[i],
        lw=5,
        xlim=(0, duration + 1),
    )
end
plt
````

save the plot
Plots.savefig(plt, "predator_prey_fig.png")

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

