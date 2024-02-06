---
title: "LinearNoiseApproximation.jl: A Julia package for linear noise approximation of stochastic reaction networks"
tags:
  - Julia
  - stochastic reaction networks
  - linear noise approximation
authors:
  - name: Xiaoming Fu
    orcid: 0000-0003-4073-9822
    affiliation: "1, 2"
  - name: Lennart Schüler
    orcid: 0000-0001-9362-1372
    affiliation: "1, 2, 3, 4"
  - name: Justin M. Calabrese
    orcid: 0000-0003-0575-6408
    affiliation: "1, 2, 5, 6"
affiliations:
  - name: Center for Advanced Systems Understanding, Görlitz, Germany
    index: 1
  - name: Helmholtz-Zentrum Dresden-Rossendorf, Dresden, Germany
    index: 2
  - name: Research Data Management - RDM, UFZ – Helmholtz Centre for Environmental Research, Leipzig, Germany
    index: 3
  - name: Dept. Monitoring and Exploration Technologies, UFZ – Helmholtz Centre for Environmental Research, Leipzig, Germany
    index: 4
  - name: Dept. of Ecological Modelling, Helmholtz Centre for Environmental Research – UFZ, Leipzig, Germany
    index: 5
  - name: Dept. of Biology, University of Maryland, College Park, MD, USA
    index: 6

date: 02 February 2024
bibliography: paper.bib
---



# Summary


Linear noise approximation (LNA) is a moment-based approximation method for stochastic chemical kinetics [@vankampenExpansionMasterEquation2007]. A general chemical kinetics reaction can be described as follows: given a set of chemical species $X_{i}, i = 1, \ldots, N$, we define $R$ chemical reactions by the notation
\begin{equation}
    \label{eq:chemical_reaction}
        k_r:\sum_{i=1}^{N}s_{ir}X_{i} \rightarrow \sum_{i=1}^{N}s'_{ir}X_{i}, \tag{1}
\end{equation}
where the stoichiometric coefficients $s_{ir}$ and $s'_{ir}$ are non-negative integer numbers denoting the numbers of reactant and product molecules, respectively. The quantity $k_r$ in \autoref{eq:chemical_reaction} is called the reaction rate constant of the $r$-th reaction. Classically, the dynamics of a chemical reaction system as in \autoref{eq:chemical_reaction} is modelled by the law of mass action. The law of mass action states that the rate of a reaction is proportional to the product of the concentrations of reactant molecules, which lead to the following rate equations as:
\begin{equation}
    \label{eq:rate_equations}
    \frac{\mathrm{d}}{\mathrm{dt}} \phi_i = \sum_{r=1}^{R} S_{ir} g_r(\boldsymbol{\phi}), \tag{2}
\end{equation}
where $\phi_i$ is the concentration of species $X_i$, 
$$
\boldsymbol{S} = \{S_{ir}\}_{N\times R},\; S_{ir}=s'_{ir} - s_{ir},\; i=1,\ldots,N,\; r=1,\ldots,R,
$$
is the stoichiometric matrix, and 
$$
g_r(\boldsymbol{\phi}) = k_r \prod_{i=1}^{N} \phi_i^{s_{ir}}, 
$$
is the rate of the $r$-th reaction. 

However, the law of mass action is only valid when the number of molecules is large. When the number of molecules is small, \autoref{eq:chemical_reaction} can instead be modelled by a continuous-time Markov jump process to study the probability of the system being in a particular state at a given time. The dynamics of such a system can be described by the chemical master equation (CME) [@schnoerrApproximationInferenceMethods2017]:
\begin{equation}
    \begin{aligned}
        \frac{\mathrm{d}}{\mathrm{dt}} P(\boldsymbol{n}, t) = \sum_{r=1}^{R} f_r(\boldsymbol{n} - \mathbf{S}_r, t) P(\boldsymbol{n} - \mathbf{S}_r, t) - \sum_{r=1}^{R} f_r(\boldsymbol{n}, t) P(\boldsymbol{n}, t), \notag
    \end{aligned}
\end{equation}
where $P(\boldsymbol{n}, t)$ is the probability of the system being in state $\boldsymbol{n}$ at time $t$, 
$$f_r(\boldsymbol{n}, t) = k_r\Omega \prod_{i=1}^{N} \frac{n_i!}{(n_i-s_{ir})!\Omega^{s_{ir}}} $$
is the propensity function of reaction $r$ at state $\boldsymbol{n}$, and $\mathbf{S}_r$ is the stoichiometric vector of reaction $r$, $\Omega$ is the system size (or volume of the system).

The chemical master equation is written directly from the rate constants and stoichiometries of all the elementary reactions of a chemical system.
When the mean and (co)variance of the system are of interest, the CME can be approximated by the LNA.
LNA is an expansion of the CME taking the inverse system size $1/\Omega$ as the perturbed variable, which is originally developed by [@**vankampenStochasticProcessesPhysics1992**]. The idea is to separate concentrations into a deterministic part, given by the solution of the deterministic rate equations, and a part describing the fluctuations about the deterministic mean
$$
\frac{n_i}{\Omega} = \phi_i  + \Omega^{-1/2} \epsilon_i,
$$
where $\phi_i$ is the solution of the deterministic rate equations (2), and $\epsilon_i$ represents fluctuations about the deterministic mean. Define $\boldsymbol{\Sigma}=\left(\epsilon_{ij}\right)_{N\times N}$ to be the covariance matrix of the fluctuations, the LNA is given by:
\begin{equation}
    \partial_t \boldsymbol{\Sigma} = \mathbf{A} \boldsymbol{\Sigma} + \boldsymbol{\Sigma} \mathbf{A}^T +  \Omega^{-1} \mathbf{B}, \tag{3}
\end{equation}
where $\mathbf{A}=\mathbf{A}(\boldsymbol{\phi}),\, \mathbf{B}=\mathbf{B}(\boldsymbol{\phi})$ both depend on the solution $\boldsymbol{\phi}$ of the rate equations \autoref{eq:rate_equations}, which are defined by 
$$ \mathbf{A} = \left(A_{ij}\right)_{N\times N}, A_{ij}=\sum_{r=1}^{R} S_{ir}\partial_{\phi_j}g_r(\boldsymbol{\phi}),
$$
as the Jacobian matrix of the deterministic rate equations, and 
$$\mathbf{B} =\left(B_{ij}\right)_{N\times N}, B_{ij} = \sum_{r=1}^{R} S_{ir}S_{jr}g_r(\boldsymbol{\phi}).$$

In this formulation, the LNA allows for analytical solutions that are locally valid close to macroscopic trajectories (solution of the rate equations) of the system. We refer to the review paper [@schnoerrApproximationInferenceMethods2017] for more details on the CME and LNA.


Deriving the LNA and applying it manually can be a cumbersome and error prone process, especially when large systems.

The solution of the CME can be computationally expensive, or even infeasible, because the set of reachable states can be huge or inﬁnite. The Linear Noise Approximation (LNA) has been introduced by Van Kampen as a second order approximation of the system size expansion of the CME (Van Kampen, 1992). It permits a stochastic characterization of the evolution of a CRN, still maintaining scalability comparable to that of the deterministic models.


# Acknowledgements

This work was partially funded by the Center of Advanced Systems Understanding (CASUS), which is financed by Germany's Federal Ministry of Education and Research (BMBF) and by the Saxon Ministry for Science, Culture and Tourism (SMWK) with tax funds on the basis of the budget approved by the Saxon State Parliament. This work was also funded by the Helmholtz Association (HGF; Helmholtz Gemeinschaft) under the framework of the "Coping capacitiy of nations facing systemic crisis - a global intercomparison exploring the SARS-CoV-2 pandemic" COCAP project with the grant number KA1-Co-10. X. F. acknowledges the support from Shanghai Sailing Program (22YF1410700).

# References

