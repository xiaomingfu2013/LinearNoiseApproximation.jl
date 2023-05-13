module LinearNoiseApproximation

using Symbolics, ModelingToolkit, Catalyst
import Catalyst: numspecies
using DiffEqBase
import DiffEqBase: ODEProblem

include("lna_expand.jl")
include("utils.jl")

export LNASystem, expand_initial_conditions
end
