"""
A type that contains the information of the LNA system.

# Fields
- `rn::ReactionSystem`: a reaction network that is used to construct the LNA system, which is from `Catalyst.jl`
- `odesys::ODESystem`: the linear noise approximation system, it is a `ODESystem` that extends the `ReactionSystem` by adding the covariance terms
- `kwargs::Any`: Optional keyword arguments
"""
struct LNASystem{Ktype}
    rn::ReactionSystem
    odesys::ODESystem
    kwargs::Ktype
end

"""
Function to obtain the ODEs according to the rate equations

# Arguments
- `rn::ReactionSystem`: a reaction network that is used to construct the LNA system, which is from `Catalyst.jl`
- `combinatoric_ratelaws::Bool`: (Optional) whether to use the combinatoric rate laws, default is `false`
"""
function get_ode_propensity(rn::ReactionSystem; combinatoric_ratelaws=false)
    return oderatelaw.(reactions(rn); combinatoric_ratelaw=combinatoric_ratelaws)
end

"""
Function to obtain the symbolic matrix of the LNA system

# Arguments
- `rn::ReactionSystem`: a reaction network that is used to construct the LNA system, which is from `Catalyst.jl`
- `combinatoric_ratelaws::Bool`: (Optional) whether to use the combinatoric rate laws, default is `false`
"""
function get_symbolic_matrix(rn::ReactionSystem; combinatoric_ratelaws=false)
    sps = Catalyst.species(rn)
    S = Catalyst.netstoichmat(rn)
    jrl = get_ode_propensity(rn; combinatoric_ratelaws=combinatoric_ratelaws)
    A_symbol = S * ModelingToolkit.jacobian(jrl, sps)
    BB_symbol = Num.(zeros(length(sps), length(sps)))
    # because the matrix is symmetric, we only need to calculate the upper triangle, and we only use the upper triangle
    N = length(sps)
    for i in 1:N, j in i:N
        BB_symbol[i, j] = sum(S[i, k] * jrl[k] * S[j, k] for k in eachindex(jrl))
    end
    return A_symbol, BB_symbol
end

"""
Function to automatically obtain the LNA system from a reaction network

# Arguments
- `rn::ReactionSystem`: a reaction network that is used to construct the LNA system, which is from `Catalyst.jl`
- `combinatoric_ratelaws::Bool`: (Optional) whether to use the combinatoric rate laws, default is `false`
"""
function LNASystem(rn::ReactionSystem; combinatoric_ratelaws=false, kwargs...)
    LNA = _get_LNA_system(rn; combinatoric_ratelaws=combinatoric_ratelaws, kwargs...)
    return LNA
end

function _get_LNA_system(rn; combinatoric_ratelaws=false, name=Symbol(string(:LNA_, rn.name)), kwargs...)
    ratesys = convert(ODESystem, rn; kwargs...)

    N = numspecies(rn)
    Σ = construct_Σ(N)
    @variables t
    A_symbol, BB_symbol = get_symbolic_matrix(rn; combinatoric_ratelaws=combinatoric_ratelaws)
    A = A_symbol * Σ
    PartialΣ = A + transpose(A) + BB_symbol
    cov_eqs = Equation[]
    for i in 1:N, j in i:N
        push!(cov_eqs, Differential(t)(Σ[i, j]) ~ PartialΣ[i, j])
    end
    connected_eqs = [equations(ratesys); cov_eqs]

    LNA = ODESystem(
        connected_eqs, t, [states(ratesys); unique(Σ)], parameters(ratesys); name=name, kwargs...
    )
    return LNASystem(rn, LNA, kwargs)
end

function Base.show(io::IO, mime::MIME"text/plain", lna::LNASystem)
    Base.show(io, mime, lna.odesys)
end

function construct_Σ(N::Int)
    @variables t
    Σ = Matrix{Num}(undef, N, N)
    for j in 1:N, i in 1:j
        Σ[i, j] = to_timedependent_symbol(:Σ, [i, j])
    end
    for j in 1:N, i in (j+1):N
        Σ[i, j] = Σ[j, i]
    end
    return Σ
end

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
    if length(u0) == numspecies(expsys.rn)
       u0_expanded = expand_initial_conditions(expsys, u0)
    elseif length(u0) == numspecies(expsys)
        u0_expanded = u0
    else
        error("the length of u0 must be equal to the number of species in the ReactionSystem or to the number of species in the LNA system")
    end
    return DiffEqBase.ODEProblem(getfield(expsys, :odesys), u0_expanded, args...; kwargs...)
end
