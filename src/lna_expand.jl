struct LNASystem{Ktype}
    Ω::Real # system size
    rn::ReactionSystem
    odesys::ODESystem
    kwargs::Ktype
end

function get_ode_propensity(rn::ReactionSystem; combinatoric_ratelaws=false)
    return oderatelaw.(reactions(rn); combinatoric_ratelaw=combinatoric_ratelaws)
end

function get_symbolic_matrix(rn::ReactionSystem; kwargs...)
    sps = Catalyst.species(rn)
    S = Catalyst.netstoichmat(rn)
    jrl = get_ode_propensity(rn; kwargs...)
    A_symbol = S * Symbolics.jacobian(jrl, sps)
    BB_symbol = Num.(zeros(length(sps), length(sps)))
    # because the matrix is symmetric, we only need to calculate the upper triangle, and we only use the upper triangle
    for j in 1:length(sps), i in 1:j
        BB_symbol[i, j] = sum(S[i, k] * jrl[k] * S[j, k] for k in eachindex(jrl))
    end
    return A_symbol, BB_symbol
end

function LNASystem(
    rn::ReactionSystem, Ω::SystemSize; combinatoric_ratelaws=false, kwargs...
) where {SystemSize<:Real}
    LNA = _get_LNA_system(rn, Ω; combinatoric_ratelaws=combinatoric_ratelaws, kwargs...)
    return LNA
end

function LNASystem(rn::ReactionSystem; combinatoric_ratelaws=false, kwargs...)
    if !(:Ω ∈ toexpr(parameters(rn)))
        @info "the system size Ω is not defined as a parameter in the Reaction system, then use the default Ω = 1"
    end
    LNA = _get_LNA_system(rn, 1; combinatoric_ratelaws=combinatoric_ratelaws, kwargs...)
    return LNA
end

function _get_LNA_system(rn, Ω; kwargs...)
    ratesys = convert(ODESystem, rn; kwargs...)

    N = numspecies(rn)
    Σ = construct_Σ(N)
    @variables t
    A_symbol, BB_symbol = get_symbolic_matrix(rn; kwargs...)
    A = A_symbol * Σ
    PartialΣ = A + transpose(A) + BB_symbol
    cov_eqs = Equation[]
    for j in 1:N, i in 1:j
        push!(cov_eqs, Differential(t)(Σ[i, j]) ~ PartialΣ[i, j])
    end
    connected_eqs = [equations(ratesys); cov_eqs]

    @named LNA = ODESystem(
        connected_eqs, t, [states(ratesys); unique(Σ)], parameters(ratesys)
    )
    return LNASystem(Ω, rn, LNA, kwargs)
end

function construct_Σ(N::Int)
    @variables t
    Σ = Matrix{Num}(undef, N, N)
    for j in 1:N, i in 1:j
        Σ[i, j] = to_timedependent_symbol(:Σ, [i, j])
    end
    for j in 1:N, i in (j + 1):N
        Σ[i, j] = Σ[j, i]
    end
    return Σ
end

"""
    function get_val_matrix_at_equilirium(rn::ReactionSystem, rates, equilibrium::Vector{T}) where {T<:Real}
"""
function get_val_matrix_at_equilirium(
    rn::ReactionSystem, rates, equilibrium::Vector{T}
) where {T<:Real}
    sps = Catalyst.species(rn)
    @assert length(sps) == length(equilibrium)
    @assert length(parameters(rn)) == length(rates)
    replace_eq = Dict(sps[i] => equilibrium[i] for i in eachindex(equilibrium))
    replace_rate = get_replacemap_rates(rn, rates)
    replace_dict = merge(replace_eq, replace_rate)
    A_symbol, BB_symbol = get_symbolic_matrix(rn)
    A = [i.val for i in substitute(A_symbol, replace_dict)]
    BB = [i.val for i in substitute(BB_symbol, replace_dict)]
    return Matrix{T}(A), Matrix{T}(BB)
end

function get_replacemap_rates(rn, rates::Vector{T}) where {T<:Real}
    return Dict(parameters(rn)[i] => rates[i] for i in eachindex(rates))
end
function get_replacemap_rates(rn, rates::Dict{Symbol,T}) where {T<:Real}
    return symmap_to_varmap(rn, rates)
end
function get_replacemap_rates(rn, rates::Vector{Pair{Num,T}}) where {T<:Real}
    return Dict(symmap_to_varmap(rn, rates))
end

function DiffEqBase.ODEProblem(expsys::LNASystem, u0, args...; kwargs...)
    @assert length(u0) == numspecies(expsys.rn) "the length of u0 must be equal to the number of species in the ReactionSystem"
    u0_expanded = expand_initial_conditions(expsys, u0)
    return DiffEqBase.ODEProblem(getfield(expsys, :odesys), u0_expanded, args...; kwargs...)
end
