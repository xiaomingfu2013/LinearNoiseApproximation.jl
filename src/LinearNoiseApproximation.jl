module LinearNoiseApproximation
using Reexport
using ModelingToolkit, Catalyst
import Catalyst: numspecies
using DiffEqBase
import DiffEqBase: ODEProblem
@reexport using Catalyst
include("lna_expand.jl")
include("utils.jl")

export LNASystem, expand_initial_conditions, find_states_cov_number
end
