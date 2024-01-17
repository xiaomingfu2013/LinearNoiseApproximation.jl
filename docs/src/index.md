```@meta
CurrentModule = LinearNoiseApproximation
```

# LinearNoiseApproximation

## Installation
```julia
using Pkg
Pkg.add("LinearNoiseApproximation")
using LinearNoiseApproximation
```
This package provides a numerical method of applying linear noise approximation (LNA) to a given reaction system using [Catalyst](https://github.com/SciML/Catalyst.jl). The derived expanded system can be solved using [DifferentialEquations.jl](https://github.com/JuliaDiffEq/DifferentialEquations.jl).

## Usage
### Reaction network with system size Ω

Suppose we have a reaction network

```julia
rn = @reaction_network TwoStage begin
    @parameters  v0 v1 d0 d1 Ω
    @species M(t) P(t)
    v0*Ω, ∅ --> M
    v1, M --> M + P
    d0, M --> ∅
    d1, P --> ∅
end
```
Here the system size `Ω` is integrated in the system (propensities), you can simply obtain the LNA system by calling
```julia
LNAsys = LNASystem(rn)
```
then solving the LNA by
```julia
@variables v0, v1, d0, d1, Ω
rates = [v0=>4.0,v1=>10.0,d0=>1.0,d1=>1.0,Ω=>5.0]
tspan = (0.0, 20.0)
u0 =[0, 0] # initial condition for M and P
prob = ODEProblem(LNAsys, u0, tspan, rates)
sol = solve(prob,Vern7(),abstol=1e-7, saveat =1.0)
```
## Brief introduction to linear noise approximation
Linear noise approximation (LNA) is a moment-based approximation method for stochastic chemical kinetics [[1]](#1). A general chemical kinetics reaction can be described as follows: given a set of chemical species $X_{i}, i = 1, \ldots, N$, we define $R$ chemical reactions by the notation
```math
\begin{equation}
        k_r:\sum_{i=1}^{N}s_{ir}X_{i} \rightarrow \sum_{i=1}^{N}s'_{ir}X_{i}, \tag{1}
\end{equation}
```
where the stoichiometric coefficients $s_{ir}$ and $s'_{ir}$ are non-negative integer numbers denoting the numbers of reactant and product molecules, respectively. The quantity $k_r$ in Equation (1) is called the reaction rate constant of the $r$-th reaction. Classically, the dynamics of a chemical reaction system as in Equation (1) is modelled by the law of mass action. The law of mass action states that the rate of a reaction is proportional to the product of the concentrations of reactant molecules, which lead to the following rate equations as:
```math
\begin{equation}
    \frac{\mathrm{d}}{\mathrm{dt}} \phi_i = \sum_{r=1}^{R} S_{ir} g_r(\boldsymbol{\phi}), \tag{2}
\end{equation}
```
where $\phi_i$ is the concentration of species $X_i$, 
```math
\boldsymbol{S} = \{S_{ir}\}_{N\times R},\; S_{ir}=s'_{ir} - s_{ir},\; i=1,\ldots,N,\; r=1,\ldots,R \notag
```
is the stoichiometric matrix, and 
```math
g_r(\boldsymbol{\phi}) = k_r \prod_{i=1}^{N} \phi_i^{s_{ir}}, 
```
is the rate of the $r$-th reaction. 

However, the law of mass action is only valid when the number of molecules is large. When the number of molecules is small, System (1) can instead be modelled by a continuous-time Markov jump process to study the probability of the system being in a particular state at a given time. The dynamics of such a system can be described by the chemical master equation (CME) [[2]](#2):
```math
\begin{equation}
    \begin{aligned}
        \frac{\mathrm{d}}{\mathrm{dt}} P(\boldsymbol{n}, t) = \sum_{r=1}^{R} f_r(\boldsymbol{n} - \mathbf{S}_r, t) P(\boldsymbol{n} - \mathbf{S}_r, t) - \sum_{r=1}^{R} f_r(\boldsymbol{n}, t) P(\boldsymbol{n}, t), \notag
    \end{aligned}
\end{equation}
```
where $P(\boldsymbol{n}, t)$ is the probability of the system being in state $\boldsymbol{n}$ at time $t$, 
```math
f_r(\boldsymbol{n}, t) = k_r\Omega \prod_{i=1}^{N} \frac{n_i!}{(n_i-s_{ir})!\Omega^{s_{ir}}} 
```
is the propensity function of reaction $r$ at state $\boldsymbol{n}$, and $\mathbf{S}_r$ is the stoichiometric vector of reaction $r$, $\Omega$ is the system size (or volume of the system).

The chemical master equation is written directly from the rate constants and stoichiometries of all the elementary reactions of a chemical system, but neither analytical nor numerical solutions are in general available. Fortunately, the chemical master equation can often be simplified in a linear noise approximation. Linear noise approximation is an expansion of the CME taking the inverse system size $1/\Omega$ as the perturbed variable, which is originally developed by [[3]](#3). The idea is to separate concentrations into a deterministic part, given by the solution of the deterministic rate equations, and a part describing the fluctuations about the deterministic mean
```math
\frac{n_i}{\Omega} = \phi_i  + \Omega^{-1/2} \epsilon_i,
```
where $\phi_i$ is the solution of the deterministic rate equations (2), and $\epsilon_i$ represents fluctuations about the deterministic mean. Define $\boldsymbol{\Sigma}=\left(\epsilon_{ij}\right)_{N\times N}$ to be the covariance matrix of the fluctuations, the linear noise approximation is given by:
```math
\begin{equation}
    \partial_t \boldsymbol{\Sigma} = \mathbf{A} \boldsymbol{\Sigma} + \boldsymbol{\Sigma} \mathbf{A}^T +  \Omega^{-1} \mathbf{B}, \tag{3}
\end{equation}
```
where $\mathbf{A}=\mathbf{A}(\boldsymbol{\phi}),\, \mathbf{B}=\mathbf{B}(\boldsymbol{\phi})$ both depend on the solution $\boldsymbol{\phi}$ of the rate equation \eqref{eq:rate_equations}, which are defined by 
```math
 \mathbf{A} = \left(A_{ij}\right)_{N\times N}, A_{ij}=\sum_{r=1}^{R} S_{ir}\partial_{\phi_j}g_r(\boldsymbol{\phi}),
```
as the Jacobian matrix of the deterministic rate equations, and 
```math
\mathbf{B} =\left(B_{ij}\right)_{N\times N}, B_{ij} = \sum_{r=1}^{R} S_{ir}S_{jr}g_r(\boldsymbol{\phi}).
```

In this formulation, the LNA allows for analytical solutions that are locally valid close to macroscopic trajectories (solution of the rate equations) of the system. We refer to the review paper [[2]](#2) for more details on the CME and LNA.

## References
<a id="1">[1]</a> Nicolaas Godfried van Kampen. The Expansion of the Master Equation. John Wiley & Sons, Inc., 2007. 

<a id="2">[2]</a> David Schnoerr, Guido Sanguinetti, and Ramon Grima. Approximation and inference methods for stochastic biochemical kinetics — a tutorial review (2017). Journal of Physics A: Mathematical and Theoretical, 50(9):093001.

<a id="3">[3]</a> Nicolaas Godfried van Kampen. Stochastic Processes in Physics and Chemistry, volume 1. Elsevier, 1992.