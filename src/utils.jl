"""
Function to expand the initial conditions of the LNA system.

# Arguments
- `sys::LNASystem`: the LNA system
- `u0::Vector`: the initial conditions of the orignal system (rate equations)
"""
function expand_initial_conditions(sys::LNASystem, u0)
    @assert length(u0) == numspecies(sys.rn) "u0 must be of the same length as the number of species"
    u0_expanded = zeros(eltype(u0), numspecies(sys))
    u0_expanded[1:length(u0)] = u0
    return u0_expanded
end

"""
Function to find the index of the covariance matrix in the LNA system.

# Arguments
- `num::Int`: the index of the species in the reaction network
- `lna_sys::LNASystem`: the LNA system
"""
function find_states_cov_number(num::Int, lna_sys::LNASystem)
    symb = to_timedependent_symbol(:Σ, [num, num])
    return num, findfirst(isequal(symb), lna_sys.odesys.states)
end

function find_states_cov_number(num::Vector{Int}, lna_sys::LNASystem)
    mean_idxs = Int[]
    var_idxs = Int[]
    for i in num
        mean_idx, var_idx = find_states_cov_number(i, lna_sys)
        push!(mean_idxs, mean_idx)
        push!(var_idxs, var_idx)
    end
    return mean_idxs, var_idxs
end


numspecies(expsys::LNASystem) = length(states(expsys.odesys))
species(expsys::LNASystem) = states(expsys.odesys)
function to_timedependent_symbol(symbol::Symbol, vec::Vector)
    @variables t
    return Num(
        ModelingToolkit.variable(symbol, join(vec); T=ModelingToolkit.FnType{Tuple{Any},Real})
    )(
        t
    )
end
