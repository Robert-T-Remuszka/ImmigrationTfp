using JLD2

include("SsSolve_Functions.jl")
include("RowCostsEstimate_Functions.jl")

# %% Structs =================================================================

"""
A perfect-foresight transition path from the pre-shock steady state (t=0) to
a horizon T assumed long enough to have returned to it, following a one-time
MIT shock to the foreign-born ROW→US migration cost. Wᵈ, Wᶠ, Lᵈ, Lᶠ, Vᵈ, Vᶠ
are N×(T+1) matrices, column k holding period t=k-1; Πᵈ, Πᶠ are N×N×(T+1).
"""
struct TransitionPath{T <: Real}
    Wᵈ::Matrix{T}
    Wᶠ::Matrix{T}
    Vᵈ::Matrix{T}
    Vᶠ::Matrix{T}
    Lᵈ::Matrix{T}
    Lᶠ::Matrix{T}
    Πᵈ::Array{T, 3}
    Πᶠ::Array{T, 3}
end

# %% Solver Functions =========================================================

"""
One backward step of the EMAX recursion: given next period's (known) value
V_next, this period's wages W and cost matrix f, compute this period's value
and choice probabilities directly. Unlike solve_value_choiceprobs, this is
not a fixed-point search -- V_next is already known (from later in the
transition, or the terminal steady state), so each location's value and
choice row is a single softmax/logsumexp evaluation.
"""
function backward_step(W::Vector{T}, f::Matrix{T}, V_next::Vector{T}, β::T, ν::T) where {T <: Real}

    N = length(W)
    lnW = log.(W)
    V = similar(V_next)
    Π = zeros(T, N, N)

    for l in 1:N
        x = (β .* V_next .- f[l, :]) ./ ν
        Π[l, :] .= softmax(x)
        V[l] = lnW[l] + ν * logsumexp(x)
    end

    return V, Π

end

"""
Foreign-born mobility-cost shock path, mₜ = ψᵗσ for t=0,...,T (a one-time
MIT shock of size σ at t=0, decaying geometrically at rate ψ thereafter --
Boppart, Krusell and Mitman (2018)'s own convention for a single AR(1)
innovation: "the shock will go up by, say, one unit in period 0 and thus
delivers the full sequence (1,ρ,ρ²,ρ³,...)"). Returned as a length-(T+1)
vector indexed the Julia way: entry k holds mₜ for t=k-1. A cost *reduction*
(this exercise) uses σ<0.
"""
mobility_cost_path(ψ::Real, σ::Real, T::Integer) = [ψ^(k - 1) * σ for k in 1:T + 1]

"""
Foreign-born migration cost matrix at a point in the transition: the
calibrated fᶠ with the mobility-cost shock mₜ added to the ROW→state margin
only (paper's ςₗₗ' indicator, eq. 2.1: row N = ROW as origin, columns 1:N-1
= US states as destination). Domestic costs and every other entry of fᶠ
(state-to-state, or state-to-ROW) are untouched by mₜ.
"""
function shocked_fᶠ(fᶠ::Matrix{T}, m_t::Real) where {T <: Real}

    N = size(fᶠ, 1)
    f = copy(fᶠ)
    f[N, 1:N - 1] .+= m_t

    return f

end

"""
Solve the perfect-foresight transition path following a one-time MIT shock
(size σ, decaying at rate ψ) to the foreign-born ROW→US migration cost,
starting from the pre-shock steady state `ss`.

Boundary conditions: labor supplies at t=0 are predetermined (workers already
committed to their t=0 location before the unanticipated shock hits, per the
model's own timing -- Lₜ₊₁ depends on choices made using period-t
information) and so equal the steady state; wages at t=0 are therefore also
unchanged from steady state, since they depend only on L₀=Lˢˢ. The shock
still affects period-0 *choices* (for period 1), since the cost term entering
that decision is fᶠ+mₜ. Wages and labor supplies are assumed to have
returned to the steady state by t=T -- a finite-horizon approximation to the
true T→∞ limit, standard in shooting algorithms; check
`Lᵈ[:,end] ≈ ss.Lᵈ` and `Lᶠ[:,end] ≈ ss.Lᶠ` on the returned path to confirm T
was large enough, and increase it if not.

Solves via damped fixed-point iteration on the wage path (guess a price path,
solve the household problem backwards from the known terminal steady state,
derive the implied aggregate/wage path forward, update the guess, repeat) --
Boppart, Krusell and Mitman (2018)'s own description of the shooting
algorithm, not a joint Newton solve, which would need a
(T-1)×102-dimensional finite-difference Jacobian, intractable at any
reasonable T.
"""
function solve_transition(p::Parameters, ss::Solution; ψ::Real, σ::Real, T::Integer,
                           damping::Real = 0.5, tol::Real = 1e-8, maxiter::Integer = 2000,
                           verbose::Bool = true)

    N = p.N
    m = mobility_cost_path(ψ, σ, T)

    # Column k holds period t=k-1. t=0 (k=1) and t=T (k=T+1) are pinned at
    # the steady state throughout -- only k=2,...,T (t=1,...,T-1) are
    # iterated on by the market-clearing update below.
    Wᵈ = repeat(ss.Wᵈ, 1, T + 1)
    Wᶠ = repeat(ss.Wᶠ, 1, T + 1)
    Vᵈ = zeros(N, T + 1)
    Vᶠ = zeros(N, T + 1)
    Πᵈ = zeros(N, N, T + 1)
    Πᶠ = zeros(N, N, T + 1)
    Lᵈ = zeros(N, T + 1)
    Lᶠ = zeros(N, T + 1)
    Lᵈ[:, 1] .= ss.Lᵈ
    Lᶠ[:, 1] .= ss.Lᶠ

    resid = Inf
    iter  = 0
    while resid > tol && iter < maxiter

        iter += 1

        # Backward pass: Vₜ known at t=T (steady state), sweep t=T-1,...,0.
        Vᵈ[:, T + 1] .= ss.Vᵈ
        Vᶠ[:, T + 1] .= ss.Vᶠ
        for k in T:-1:1
            Vᵈ[:, k], Πᵈ[:, :, k] = backward_step(Wᵈ[:, k], p.fᵈ, Vᵈ[:, k + 1], p.β, p.νᵈ)
            Vᶠ[:, k], Πᶠ[:, :, k] = backward_step(Wᶠ[:, k], shocked_fᶠ(p.fᶠ, m[k]), Vᶠ[:, k + 1], p.β, p.νᶠ)
        end

        # Forward pass: L₀=Lˢˢ known, sweep t=1,...,T using the Πₜ's just computed.
        for k in 1:T
            Lᵈ[:, k + 1] = Πᵈ[:, :, k]' * Lᵈ[:, k]
            Lᶠ[:, k + 1] = Πᶠ[:, :, k]' * Lᶠ[:, k]
        end

        # Market clearing at each interior period (k=2,...,T; the endpoints
        # are pinned, not updated) gives the implied wage path; ROW's wage
        # (row N) stays fixed throughout, exactly as in the steady state.
        Wᵈ_new = copy(Wᵈ)
        Wᶠ_new = copy(Wᶠ)
        for k in 2:T, l in 1:N - 1
            Wᵈ_new[l, k], Wᶠ_new[l, k] = market_clearing_wage(Lᵈ[l, k], Lᶠ[l, k], p)
        end

        resid = max(maximum(abs, log.(Wᵈ_new[:, 2:T]) .- log.(Wᵈ[:, 2:T])),
                    maximum(abs, log.(Wᶠ_new[:, 2:T]) .- log.(Wᶠ[:, 2:T])))

        Wᵈ[:, 2:T] .= damping .* Wᵈ_new[:, 2:T] .+ (1 - damping) .* Wᵈ[:, 2:T]
        Wᶠ[:, 2:T] .= damping .* Wᶠ_new[:, 2:T] .+ (1 - damping) .* Wᶠ[:, 2:T]

        verbose && iter % 25 == 0 && (println("shooting iter $iter: max log-wage resid = $resid"); flush(stdout))

    end

    resid > tol && error("Shooting algorithm failed to converge in $maxiter iterations (resid=$resid)")
    verbose && println("Shooting converged in $iter iterations, resid=$resid")

    return TransitionPath(Wᵈ, Wᶠ, Vᵈ, Vᶠ, Lᵈ, Lᶠ, Πᵈ, Πᶠ)

end

# %% Helper Functions =========================================================

"""
National aggregates -- Z, the task-weighted labor aggregate L, and
employment-weighted average wages by nativity -- from a cross-section of US
state-level wages and labor supplies. These are the four objects MakeIRF.do's
empirical LPIV estimates (Z, L, Wage_Foreign, Wage_Domestic), computed the
same way from the steady state and from every period of a transition path so
they're directly comparable. Only states 1:N-1 enter (Rest of World, index
N, is excluded, exactly as in the empirical state panel). Z and wages are
employment-weighted averages across states, matching the `wt(emp)` weighting
in the empirical `ivreg2` regression; L is a national total, since it's
already additive across states by construction (eq. 2.6).
"""
function national_aggregates(Wᵈ::Vector{T}, Wᶠ::Vector{T}, Lᵈ::Vector{T}, Lᶠ::Vector{T}, p::Parameters) where {T <: Real}

    Nstates = p.N - 1
    Z = zeros(T, Nstates)
    L = zeros(T, Nstates)
    for l in 1:Nstates
        w = Wᵈ[l] / Wᶠ[l]
        agg = TaskAggregates_LN(p.ρ, p.μ_z, p.ξ_ω, p.ξ_z, w)
        Z[l] = agg.Z
        L[l] = LaborAggregate(agg.λ, p.ρ, Lᶠ[l], Lᵈ[l], agg.oneMinusλ)
    end

    Z_nat  = sum(Z .* L) / sum(L)
    L_nat  = sum(L)
    WF_nat = sum(Wᶠ[1:Nstates] .* Lᶠ[1:Nstates]) / sum(Lᶠ[1:Nstates])
    WD_nat = sum(Wᵈ[1:Nstates] .* Lᵈ[1:Nstates]) / sum(Lᵈ[1:Nstates])

    return (; Z_nat, L_nat, WF_nat, WD_nat)

end

"""
National-aggregate IRFs (log deviation from steady state) for Z, L,
Wage_Foreign, Wage_Domestic at every horizon of a solved transition path --
see national_aggregates for how each is constructed. Returns a NamedTuple of
four length-(T+1) vectors, matching MakeIRF.do's `Z Wage_Domestic
Wage_Foreign L` outcome set and combined-graph order.
"""
function national_irfs(path::TransitionPath, ss::Solution, p::Parameters)

    T1 = size(path.Wᵈ, 2)
    ss_agg = national_aggregates(ss.Wᵈ, ss.Wᶠ, ss.Lᵈ, ss.Lᶠ, p)

    Z  = zeros(T1); L = zeros(T1); WF = zeros(T1); WD = zeros(T1)
    for k in 1:T1
        agg = national_aggregates(path.Wᵈ[:, k], path.Wᶠ[:, k], path.Lᵈ[:, k], path.Lᶠ[:, k], p)
        Z[k]  = log(agg.Z_nat  / ss_agg.Z_nat)
        L[k]  = log(agg.L_nat  / ss_agg.L_nat)
        WF[k] = log(agg.WF_nat / ss_agg.WF_nat)
        WD[k] = log(agg.WD_nat / ss_agg.WD_nat)
    end

    return (; Z, L, WF, WD)

end
