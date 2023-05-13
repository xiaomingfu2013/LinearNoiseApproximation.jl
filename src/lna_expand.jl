struct LNASystem{Ktype}
    Ω::Real # system size
    rn::ReactionSystem
    odesys::ODESystem
    kwargs::Ktype
end

get_ode_propensity(rn::ReactionSystem;combinatoric_ratelaws = false) = oderatelaw.(reactions(rn); combinatoric_ratelaw = combinatoric_ratelaws)

function get_symbolic_matrix(rn::ReactionSystem; kwargs...)
    sps = Catalyst.species(rn)
    S = Catalyst.netstoichmat(rn)
    jrl = get_ode_propensity(rn; kwargs...)
    A_symbol = S*Symbolics.jacobian(jrl, sps)
    f_diag_symbol = Num.(zeros(length(jrl),length(jrl)))
    [f_diag_symbol[i,i] = jrl[i] for i in eachindex(jrl)]
    BB_symbol = S*f_diag_symbol*transpose(S)
    return A_symbol, BB_symbol
end


function LNASystem(rn::ReactionSystem, Ω::SystemSize; combinatoric_ratelaws = false, kwargs...) where {SystemSize<:Real}
    LNA = _get_LNA_system(rn, Ω; combinatoric_ratelaws = combinatoric_ratelaws, kwargs...)
    return LNA
end

function LNASystem(rn::ReactionSystem; combinatoric_ratelaws = false, kwargs...)
    if !(:Ω ∈ toexpr(parameters(rn)))
        @warn "the system size Ω is not defined as a parameter in the Reaction system, then use the default Ω = 1"
    end
    LNA = _get_LNA_system(rn, 1.; combinatoric_ratelaws = combinatoric_ratelaws, kwargs...)
    return LNA
end

function _get_LNA_system(rn, Ω;  kwargs...)
    ratesys = convert(ODESystem, rn; kwargs...)

    N = numspecies(rn)
    Σ = construct_Σ(N)
    @variables t
    A_symbol, BB_symbol = get_symbolic_matrix(rn; kwargs...)

    PartialΣ = A_symbol*Σ+Σ*transpose(A_symbol) + BB_symbol*Ω
    cov_eqs = vcat([[Differential(t)(Σ[i,j]) ~ PartialΣ[i,j] for i in 1:j] for j in 1:N]...)

    connected_eqs = [equations(ratesys);cov_eqs]

    @named LNA = ODESystem(connected_eqs,t,[states(ratesys);unique(Σ)],parameters(ratesys))
    LNASystem(Ω, rn, LNA, kwargs)
end

function construct_Σ(N::Int)
    @variables t
    Σ = Matrix{Any}(undef, N, N)
    for j in 1:N, i in 1:j
        Σ[i,j] =  to_timedependent_symbol(:Σ, [i, j])
    end
    for j in 1:N, i in j+1:N
        Σ[i,j] = Σ[j,i]
    end
    return Σ
end



"""
    function get_val_matrix_at_equilirium(rn::ReactionSystem, rates, equilibrium::Vector{T}) where {T<:Real}
"""
function get_val_matrix_at_equilirium(rn::ReactionSystem, rates, equilibrium::Vector{T2}) where {T2<:Real}
    sps = Catalyst.species(rn)
    @assert length(sps) == length(equilibrium)
    @assert length(parameters(rn)) == length(rates)
    replace_eq = Dict(sps[i]=>equilibrium[i] for i in eachindex(equilibrium))
    replace_rate = get_replacemap_rates(rn, rates)
    T1 = valtype(replace_rate)
    T = typeintersect(T1, T2)
    replace = merge(replace_eq, replace_rate)
    A_symbol, BB_symbol = get_symbolic_matrix(rn)
    A = [i.val  for i in substitute(A_symbol, replace)]
    BB = [i.val for i in substitute(BB_symbol, replace)]
    return Matrix{T}(A), Matrix{T}(BB)
end


function get_replacemap_rates(rn, rates::Vector{T}) where {T<:Real}
    Dict(parameters(rn)[i]=>rates[i] for i in eachindex(rates))
end
function get_replacemap_rates(rn, rates::Dict{Symbol, T}) where {T<:Real}
    symmap_to_varmap(rn, rates)
end
function get_replacemap_rates(rn, rates::Vector{Pair{Num,T}}) where {T<:Real}
    Dict(symmap_to_varmap(rn, rates))
end

function DiffEqBase.ODEProblem(expsys::LNASystem, args...; kwargs...)
    DiffEqBase.ODEProblem(getfield(expsys,:odesys), args...; kwargs...)
end
