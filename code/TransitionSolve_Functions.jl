using JLD2

include("SsSolve_Functions.jl")
include("RowCostsEstimate_Functions.jl")

# %% Structs =================================================================

"""
A perfect-foresight transition path from a steady state (t=0) to a horizon T,
following an unanticipated change in the foreign-born migration cost path.
Wᵈ, Wᶠ, Lᵈ, Lᶠ, Vᵈ, Vᶠ
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
Solve a perfect-foresight transition path from the steady state `ss`, given a
path of foreign-born migration costs: `fᶠ_path` is N×N×(T+1), with slice k
holding the cost matrix of period t=k-1. The path is unanticipated at t=0.

Boundary conditions: labor supplies at t=0 are predetermined (workers already
committed to their t=0 location before the news arrives, per the model's own
timing -- Lₜ₊₁ depends on choices made using period-t information) and so
equal the steady state; wages at t=0 are therefore also unchanged from steady
state, since they depend only on L₀=Lˢˢ. The new costs still affect period-0
*choices* (for period 1). Terminal values at t=T hold the market-clearing
wages at L_T and the costs at t=T constant forever, so no terminal steady
state is needed as an input. This also covers the closed border counterfactual, 
where the foreign wage in a US state rises without bound as Lᶠ → 0.

Solves via damped fixed-point iteration on the wage path (guess a price path,
solve the household problem backwards from the terminal values, derive the
implied aggregate/wage path forward, update the guess, repeat).
"""
function solve_transition(p::Parameters, ss::Solution, fᶠ_path::Array{<:Real, 3};
                           damping::Real = 0.5, tol::Real = 1e-8, maxiter::Integer = 2000,
                           verbose::Bool = true)

    N = p.N
    T = size(fᶠ_path, 3) - 1

    # Column k holds period t=k-1. t=0 (k=1) is pinned at the steady state
    # throughout -- only k=2,...,T+1 (t=1,...,T) are iterated on by the
    # market-clearing update below.
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

        # Backward pass: terminal values hold t=T wages and costs constant,
        # then sweep t=T-1,...,0.
        Vᵈ[:, T + 1], Πᵈ[:, :, T + 1] = solve_value_choiceprobs(Wᵈ[:, T + 1], p.fᵈ, p.β, p.νᵈ)
        Vᶠ[:, T + 1], Πᶠ[:, :, T + 1] = solve_value_choiceprobs(Wᶠ[:, T + 1], fᶠ_path[:, :, T + 1], p.β, p.νᶠ)
        for k in T:-1:1
            Vᵈ[:, k], Πᵈ[:, :, k] = backward_step(Wᵈ[:, k], p.fᵈ, Vᵈ[:, k + 1], p.β, p.νᵈ)
            Vᶠ[:, k], Πᶠ[:, :, k] = backward_step(Wᶠ[:, k], fᶠ_path[:, :, k], Vᶠ[:, k + 1], p.β, p.νᶠ)
        end

        # Forward pass: L₀=Lˢˢ known, sweep t=1,...,T using the Πₜ's just computed.
        for k in 1:T
            Lᵈ[:, k + 1] = Πᵈ[:, :, k]' * Lᵈ[:, k]
            Lᶠ[:, k + 1] = Πᶠ[:, :, k]' * Lᶠ[:, k]
        end

        # Market clearing at each period after t=0 gives the implied wage path;
        # ROW's wage (row N) stays fixed throughout, exactly as in the steady state.
        Wᵈ_new = copy(Wᵈ)
        Wᶠ_new = copy(Wᶠ)
        for k in 2:T + 1, l in 1:N - 1
            Wᵈ_new[l, k], Wᶠ_new[l, k] = market_clearing_wage(Lᵈ[l, k], Lᶠ[l, k], p)
        end

        resid = max(maximum(abs, log.(Wᵈ_new[:, 2:T + 1]) .- log.(Wᵈ[:, 2:T + 1])),
                    maximum(abs, log.(Wᶠ_new[:, 2:T + 1]) .- log.(Wᶠ[:, 2:T + 1])))

        Wᵈ[:, 2:T + 1] .= damping .* Wᵈ_new[:, 2:T + 1] .+ (1 - damping) .* Wᵈ[:, 2:T + 1]
        Wᶠ[:, 2:T + 1] .= damping .* Wᶠ_new[:, 2:T + 1] .+ (1 - damping) .* Wᶠ[:, 2:T + 1]

        verbose && iter % 25 == 0 && (println("shooting iter $iter: max log-wage resid = $resid"); flush(stdout))

    end

    resid > tol && error("Shooting algorithm failed to converge in $maxiter iterations (resid=$resid)")
    verbose && println("Shooting converged in $iter iterations, resid=$resid")

    return TransitionPath(Wᵈ, Wᶠ, Vᵈ, Vᶠ, Lᵈ, Lᶠ, Πᵈ, Πᶠ)

end
