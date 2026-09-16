#================================================================
SMM simulator building blocks: BKM (2018) impulse-response-as-numerical-derivative
kernel construction around the model's own steady state.

Requires Globals.jl, Estimation_Funcs.jl, AggSupply_Functions.jl,
Solve_Baseline_Functions.jl, Solve_Counterfactual_Functions.jl to already be included.
================================================================#

#================================================================
                        STEADY STATE
================================================================#
"""
Per-period max proportional deviation from flat (period t+1 vs period t) for wages and
labor stocks, t = 1,…,T-1. A genuine steady state requires all four to be ~0 for every
period from some t★ onward -- Π's own flatness follows from L's (L_{t+1}=Πₜ'Lₜ) and from
SolveBaseline's own U̇→1 convergence criterion, so it is not tracked separately here.
"""
function level_convergence(S::Soln)

    dev(X) = [maximum(abs.(X[:, t + 1] ./ X[:, t] .- 1)) for t in 1:size(X, 2) - 1]

    return (; Wᵈ = dev(S.Wᵈ), Wᶠ = dev(S.Wᶠ), Lᵈ = dev(S.Lᵈ), Lᶠ = dev(S.Lᶠ))

end

"""
Smallest period t★ such that every series in level_convergence stays below tol for all
t ∈ [t★, T-1] -- i.e. the solution has genuinely settled, not merely touched tol once.
Returns nothing if no such t★ exists.
"""
function steady_state_period(S::Soln; tol::Real = 1e-8)

    lc    = level_convergence(S)
    worst = [max(lc.Wᵈ[t], lc.Wᶠ[t], lc.Lᵈ[t], lc.Lᶠ[t]) for t in 1:S.T - 1]

    for t★ in 1:S.T - 1
        all(worst[t★:end] .< tol) && return t★
    end

    return nothing

end

"""
Max proportional deviation from flat anywhere in S -- the scalar assertion used to certify
a re-anchored Baseline_ss is genuinely a fixed point.
"""
function flatness(S::Soln)

    lc = level_convergence(S)
    return max(maximum(lc.Wᵈ), maximum(lc.Wᶠ), maximum(lc.Lᵈ), maximum(lc.Lᶠ))

end

"""
Build a new Parameters pair anchored at period t of S instead of 1996: the labor stocks and
lagged choice probabilities at t become the new Lᵈ₀/Lᶠ₀/Πᵈ₋/Πᶠ₋, and Y₀ is backed out by
inverting resource feasibility, (1-θ)Y=wᵈlᵈ+wᶠlᶠ, at S's period-t wages/stocks -- this makes
the existing, unmodified solve_initial_wages reproduce those same wages exactly, so no solver
code changes are needed to re-anchor the model at its own steady state.
"""
function steady_state_parameters(S::Soln, t::Integer; p::Parameters)

    (; N, θ) = p
    Wᵈ_t, Wᶠ_t = S.Wᵈ[:, t], S.Wᶠ[:, t]
    Lᵈ_t, Lᶠ_t = S.Lᵈ[:, t], S.Lᶠ[:, t]
    Πᵈ_lag = t == 1 ? p.Πᵈ₋ : S.Πᵈ[:, :, t - 1]
    Πᶠ_lag = t == 1 ? p.Πᶠ₋ : S.Πᶠ[:, :, t - 1]

    Y₀ = fill(NaN, N)
    for l in 1:N - 1
        Y₀[l] = (Wᵈ_t[l] * Lᵈ_t[l] + Wᶠ_t[l] * Lᶠ_t[l]) / (1 - θ[l])
    end

    p_ss = Parameters(p.β, p.r, p.δ, p.ρ, p.θ, p.μ_z, p.ξ_ω, p.ξ_z, p.ψ, p.νᵈ, p.νᶠ,
                       N, p.wᵈ_row, p.wᶠ_row, Πᵈ_lag, Πᶠ_lag, Lᵈ_t, Lᶠ_t, Y₀)

    return p_ss, Wᵈ_t, Wᶠ_t

end

"""
Solve the model out to its own long-run steady state, then re-anchor a fresh Parameters/Soln
pair AT that steady state rather than at the 1996 transition. BKM's method linearizes around
a steady state; shocking the raw 1996-anchored transition instead gives a kernel whose shape
depends on when along the transition the shock hits (it fails time-invariance), and firing
the shock only partway through a *solved* transition would additionally let perfect-foresight
agents see it coming and pre-migrate, contaminating even the k=0 response. Re-anchoring avoids
both: take the transition's terminal state, once it has genuinely gone flat (not merely
SolveBaseline's own U̇→1 criterion -- see steady_state_period), and rebuild Parameters with
that state as its own "period zero" (steady_state_parameters). Refines by re-anchoring on the
re-solved Baseline's own terminal state, since re-solving from an already-near-fixed-point
guess is far cheaper than extending the original transition indefinitely.
"""
function SolveSteadyState(p::Parameters;
    T0::Integer = 25, T_step::Integer = 50, level_tol::Real = 1e-4, max_T::Integer = 1500,
    T_ss::Integer = 120, refine_maxiter::Integer = 10,
    outer_tol::Real = 1e-6, outer_maxiter::Integer = 20_000, ss_tol::Real = 1e-6,
    verbose::Bool = true)

    Transition = SolveBaseline(p; T0, outer_tol, outer_maxiter, ss_tol, verbose)
    t★ = steady_state_period(Transition; tol = level_tol)

    while isnothing(t★) && Transition.T < max_T
        verbose && println("Extending T from $(Transition.T) to $(Transition.T + T_step) to reach level_tol=$level_tol")
        Transition = SolveBaseline(p; init = ExtendSoln(Transition, Transition.T + T_step),
                                    outer_tol, outer_maxiter, ss_tol, verbose)
        t★ = steady_state_period(Transition; tol = level_tol)
    end
    isnothing(t★) && error("Transition did not reach level_tol=$level_tol within max_T=$max_T periods")

    p_ss, _, _  = steady_state_parameters(Transition, t★; p)
    Baseline_ss = SolveBaseline(p_ss; T0 = T_ss, outer_tol, outer_maxiter, ss_tol, verbose)

    iters = 0
    while flatness(Baseline_ss) >= level_tol && iters < refine_maxiter
        p_ss, _, _  = steady_state_parameters(Baseline_ss, Baseline_ss.T; p = p_ss)
        Baseline_ss = SolveBaseline(p_ss; T0 = T_ss, outer_tol, outer_maxiter, ss_tol, verbose)
        iters += 1
    end

    fl = flatness(Baseline_ss)
    fl < level_tol || error("SolveSteadyState failed to reach flatness < $level_tol after $refine_maxiter refinements (flatness = $fl)")

    return (; p_ss, Baseline_ss, Transition, t★, flatness = fl, iters)

end

#================================================================
                        BKM KERNEL
================================================================#
"""
BKM (2018) impulse-response-as-numerical-derivative kernel: solve one perfect-foresight
counterfactual response to a probe-sized MIT shock (mit_shock(σ = σ_p, ψ)) fired at t=1
against the (already re-anchored, flat) steady-state Baseline, then normalize by σ_p to get
the "response per unit shock" 𝒦(k) = ln(x^CF/x^Baseline)[k] / σ_p for k = 0,…,K-1 (k=0 is the
first period the shock bites, Soln column 2), for the 51 US locations (excludes Rest-of-World).

σ_p > 0 is a mobility-COST INCREASE (cost_matrix multiplies the ROW→US entries, and
probabilities enter as C^(-1/ν)), so fewer ROW→US migrants -- confirmed against
CounterfactualIRF.pdf. σ_p is a purely numerical probe size (a finite-difference step), not
the structural innovation SD σ -- see kernel_scalability for the check that a given σ_p is
small enough to trust as a linear approximation.
"""
function bkm_kernel(Baseline::Soln; p::Parameters, σ_p::Real, K::Integer,
    CF_tol::Real = 1e-6, CF_maxiter::Integer = 20_000, verbose::Bool = false)

    @assert Baseline.T >= K + 1 "Baseline horizon T=$(Baseline.T) too short for K=$K"

    ψ = p.ψ
    M̂ = [mit_shock(; σ = σ_p, ψ)(t) for t in 1:Baseline.T - 1]
    CF = SolveCounterfactual(M̂; Baseline, p, CF_tol, CF_maxiter, verbose)

    us = 1:p.N - 1
    Z_B, Z_CF = ComputeZ(Baseline; p), ComputeZ(CF; p)
    L_B, L_CF = ComputeL(Baseline; p), ComputeL(CF; p)

    κ = (
        Z  = log.(Z_CF[us, 2:K + 1]  ./ Z_B[us, 2:K + 1])  ./ σ_p,
        L  = log.(L_CF[us, 2:K + 1]  ./ L_B[us, 2:K + 1])  ./ σ_p,
        Wᵈ = log.(CF.Wᵈ[us, 2:K + 1] ./ Baseline.Wᵈ[us, 2:K + 1]) ./ σ_p,
        Wᶠ = log.(CF.Wᶠ[us, 2:K + 1] ./ Baseline.Wᶠ[us, 2:K + 1]) ./ σ_p,
    )

    return κ, CF

end

#================================================================
                        DIAGNOSTICS
================================================================#
"""
BKM's own accuracy check (their §5.2, Figs. 5-6): solve the kernel at several probe sizes
(and signs) and compare the *normalized* kernels 𝒦(σ_p) = IRF(σ_p)/σ_p against a small
reference probe. Under exact linearity 𝒦 would not depend on σ_p at all; the report's
`rel_dist` column is how far that fails at each probe size, and `antisym_dist` checks
𝒦(σ_p) ≈ -𝒦(-σ_p) (mirror-image responses under linearity). Pick the largest σ_p whose
rel_dist stays within tolerance -- that is the probe with the best signal-to-noise that is
still safely inside the linear region.
"""
function kernel_scalability(Baseline::Soln; p::Parameters, K::Integer,
    probes::Vector{<:Real} = [0.01, 0.05, 0.10, 0.25, 0.50, 1.00, -0.01, -0.05, -0.10, -0.25, -0.50, -1.00],
    ref::Real = 0.01, CF_tol::Real = 1e-6, CF_maxiter::Integer = 20_000, verbose::Bool = false)

    @assert ref in probes "ref=$ref must be included in probes"

    kernels = Dict(σ_p => first(bkm_kernel(Baseline; p, σ_p, K, CF_tol, CF_maxiter, verbose)) for σ_p in probes)
    κ_ref = kernels[ref]

    rows = NamedTuple[]
    for σ_p in probes, outcome in (:Z, :L, :Wᵈ, :Wᶠ)

        d   = maximum(abs.(kernels[σ_p][outcome] .- κ_ref[outcome]))
        rel = d / maximum(abs.(κ_ref[outcome]))

        antisym = haskey(kernels, -σ_p) ?
            maximum(abs.(kernels[σ_p][outcome] .+ kernels[-σ_p][outcome])) : missing

        push!(rows, (; σ_p, outcome, sup_dist = d, rel_dist = rel, antisym_dist = antisym))

    end

    return kernels, DataFrame(rows)

end
