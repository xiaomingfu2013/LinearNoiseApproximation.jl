function expand_initial_conditions(sys::LNASystem, u0)
    @assert length(u0) == numspecies(sys.rn) "u0 must be of the same length as the number of species"
    u0_expanded = zeros(eltype(u0), numspecies(sys))
    u0_expanded[1:length(u0)] = u0
    return u0_expanded
end
numspecies(expsys::LNASystem) = length(states(expsys.odesys))
species(expsys::LNASystem) = states(expsys.odesys)
to_timedependent_symbol(symbol::Symbol, vec::Vector) =(@variables t; Num(Symbolics.variable(symbol, join(vec); T = ModelingToolkit.FnType{Tuple{Any},Real}))(t))

function find_states_cov_number(num::Int, lna_sys::LNASystem)
    symb = to_timedependent_symbol(:Σ, [num, num])
    return num, findfirst(isequal(symb), lna_sys.odesys.states)
end
