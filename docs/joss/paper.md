---
title: 'LinearNoiseApproximation.jl: A Julia package for linear noise approximation of stochastic reaction networks'
tags:
  - Julia
  - reaction-systems
  - linear-noise-approximation
authors:
  - name: Xiaoming Fu
    orcid: 0000-0003-4073-9822
    affiliation: "1, 2"
  - name: Lennart Schüler
    orcid: 0000-0001-9362-1372
    affiliation: "3, 4"
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
date: 08 February 2024
bibliography: paper.bib
---
# Summary

Stochastic reaction networks (SRNs) were originally developed in the 1960s to
model chemical systems based on the law of mass action. However, SRNs have
proven their flexibility beyond chemical systems and across many different
fields. These include biochemistry with gene expression models, epidemiology
with compartmental models, ecology with the Lotka-Volterra equations,
condensed matter and semiconductor physics with diffusion processes, and even
sociology [@voit_150_2015; @saa_formulation_2017; @yu_mathematical_2018;
@feinberg_foundations_2019; @gorban_quasichemical_2011;
@galam_sociophysics_2012] to name just a few.

SRNs include stochastic fluctuations and are often non-linear. Despite their
broad utility, direct numerical computation of SRNs often becomes intractable
for large system sizes, which can already happen with only a few "components",
i.e. reactants, genes, compartments, ... [@schnoerr_approximation_2017].
However, in many applications of SRNs, we are only interested in the first
moments of the system, e.g. the mean and the (co)variance, which can be
approximated and efficiently calculated with the widely-used Linear Noise
Approximation (LNA), derived through a series expansion. Here, we introduce
LinearNoiseApproximation.jl, which is a Julia package that automates LNA
derivation in a simple, robust, and general way and integrates seamlessly
with popular scientific software packages available for Julia.

# Statement of need

As described above, the LNA is very widely used across the sciences. Deriving
the LNA manually, however, is tedious and error prone [@fu_simultaneous_2024].
Automated approaches to LNA derviation have therefore been developed
[@fan_means_2016; @thomas_intrinsic_2012], but neither of these existing
packages is actively maintained, and both have significant usability
limitations. Specifically, the Python package MEANS [@fan_means_2016] has not
been updated since 2016 and relies on Python support for solving differential
equations, which is arguably less robust compared to the differential equation
solvers available for Julia, especially for stochastic differential equations
[@rackauckas_differential_2017]. Similarly, the C++ package iNA
[@thomas_intrinsic_2012] was last updated 2014 and has to be used via its GUI
or C++ directly, making it less compatible with modern scientific software
stacks. To the best of our knowledge, no other open-source software packages
for automated LNA derivation exist, which represents a significant tool gap
for researchers who rely on SRNs to model dynamical systems across a broad
array of scientific disciplines.

We believe that Julia is a natural choice as a programming language for
automated LNA derivation. The Julia ecosystem includes Catalyst.jl
[@loman_catalyst_2023], making it easy and intuitive to construct a reaction
network, which our package can then automatically translate to a linear noise
approximated ODE system on-the-fly. This system in turn can then be solved
with the Julia package DifferentialEquations.jl
[@rackauckas_differential_2017].
@fu_simultaneous_2024 use LinearNoiseApproximation.jl to simultaneously
estimate model parameter and changepoints at which model parameters change in
switching dynamical systems. This is done over a broad range of
discipline-specific examples including a highly oscillatory non-linear
predator-prey model, a stochastic gene expression model, and a
compartment-type epidemiological model with only limited and noisy data
available.

# Mathematics

While LNA has applications in many different fields, it was first derived as a
moment-based approximation method for stochastic chemical kinetics
[@van_kampen_expansion_2007]. For brevity and historical consistency, we will
therefore stick to the naming conventions from this field. Chemical kinetics
reactions with a small number of molecules can be modeled by a continuous-time
Markov jump process to study the probability of the system being in a
particular state at a given time. The dynamics of such a system can be
described by the chemical master equation (CME) [@schnoerr_approximation_2017]:

\begin{equation}
  \begin{aligned}
    \frac{\mathrm{d}}{\mathrm{dt}} P(\boldsymbol{n}, t) =
    \sum_{r=1}^{R} f_r(\boldsymbol{n} - \mathbf{S}_r, t) P(\boldsymbol{n} -
    \mathbf{S}_r, t) - \sum_{r=1}^{R} f_r(\boldsymbol{n}, t)
    P(\boldsymbol{n}, t) \; ,
  \end{aligned}
\end{equation}

where $P(\boldsymbol{n}, t)$ is the probability of the system being in state
$\boldsymbol{n}$ at time $t$. $f_r$ is the propensity function of reaction
$r$ at state $\boldsymbol{n}$, $\mathbf{S}_r$ is the stoichiometric vector
of reaction $r$.
The CME is written directly from the rate constants and stoichiometries of all
the elementary reactions of a chemical system. When the mean and (co)variance
of the system are of interest, the CME can be approximated by the LNA. LNA is
an expansion of the CME taking the inverse system size $1/\Omega$ as the
perturbed variable, which was originally developed by
[@van_kampen_stochastic_1992]. The idea is to separate concentrations into a
deterministic part, given by the solution of the deterministic rate equations,
and a part describing the fluctuations about the deterministic mean

\begin{equation}
  \frac{n_i}{\Omega} = \phi_i  + \Omega^{-1/2} \epsilon_i,
\end{equation}

where $\phi_i$ is the concentration of species $i$, and $\epsilon_i$
represents fluctuations about the deterministic mean. Define
$\boldsymbol{\Sigma}=\left(\epsilon_{ij}\right)_{N\times N}$ to be the
covariance matrix of the fluctuations and the LNA is given by:

\begin{equation}
  \partial_t \boldsymbol{\Sigma} = \mathbf{A} \boldsymbol{\Sigma} +
  \boldsymbol{\Sigma} \mathbf{A}^T +  \Omega^{-1} \mathbf{B} \; ,
\end{equation}

with $\mathbf{A}(\boldsymbol{\phi})$ and $\mathbf{B}(\boldsymbol{\phi})$ being
defined by 

\begin{equation}
  \mathbf{A} = \left(A_{ij}\right)_{N\times N}, A_{ij}=\sum_{r=1}^{R}
  S_{ir}\partial_{\phi_j}g_r(\boldsymbol{\phi}) \; ,
\end{equation}

as the Jacobian matrix of the deterministic rate equations and 

\begin{equation}
  \mathbf{B} =\left(B_{ij}\right)_{N\times N}, B_{ij} =
  \sum_{r=1}^{R} S_{ir}S_{jr}g_r(\boldsymbol{\phi}).
\end{equation}

In this formulation, the LNA allows for analytical solutions that are locally
valid close to macroscopic trajectories (solution of the rate equations) of
the system. We refer to the review paper [@schnoerr_approximation_2017]
for more details on the CME and LNA.

# Acknowledgements

This work was partially funded by the Center of Advanced Systems Understanding
(CASUS), which is financed by Germany's Federal Ministry of Education and
Research (BMBF) and by the Saxon Ministry for Science, Culture and Tourism
(SMWK) with tax funds on the basis of the budget approved by the Saxon State
Parliament. This work was also funded by the Helmholtz Association
(HGF; Helmholtz Gemeinschaft) under the framework of the "Coping capacity of
nations facing systemic crisis - a global intercomparison exploring the
SARS-CoV-2 pandemic" COCAP project with the grant number KA1-Co-10. X.F.
acknowledges the support from Shanghai Sailing Program (22YF1410700).

# References
