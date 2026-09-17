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
    T_ss::Integer = 400, refine_maxiter::Integer = 10,
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

Returns kernels for the PRIMITIVES (Lᵈ, Lᶠ, Wᵈ, Wᶠ) -- these are what simulate_panel should
convolve against drawn innovations -- plus Z and L (the CES task-productivity/labor
aggregates) computed by directly log-linearizing ComputeZ/ComputeL's output, kept only as a
diagnostic against construct_aggregates's nonlinear reconstruction from the simulated
primitives (see its docstring for why the two are not interchangeable for building a
simulated panel).

σ_p > 0 is a mobility-COST INCREASE (cost_matrix multiplies the ROW→US entries, and
probabilities enter as C^(-1/ν)), so fewer ROW→US migrants -- confirmed against
CounterfactualIRF.pdf. σ_p is a purely numerical probe size (a finite-difference step), not
the structural innovation SD σ -- see kernel_scalability for the check that a given σ_p is
small enough to trust as a linear approximation.

The ψ keyword overrides p.ψ without rebuilding Parameters -- safe because ψ is read in
exactly one place in the whole solver (mit_shock); Baseline itself never depends on ψ
(SolveBaseline/SolveSteadyState never read p.ψ), so the SAME cached Baseline is valid for
every ψ in a (ψ,σ) grid search -- only this kernel call needs rebuilding per ψ.
"""
function bkm_kernel(Baseline::Soln; p::Parameters, σ_p::Real, K::Integer, ψ::Real = p.ψ,
    CF_tol::Real = 1e-6, CF_maxiter::Integer = 20_000, verbose::Bool = false)

    @assert Baseline.T >= K + 1 "Baseline horizon T=$(Baseline.T) too short for K=$K"

    M̂ = [mit_shock(; σ = σ_p, ψ)(t) for t in 1:Baseline.T - 1]
    CF = SolveCounterfactual(M̂; Baseline, p, CF_tol, CF_maxiter, verbose)

    us = 1:p.N - 1
    Z_B, Z_CF = ComputeZ(Baseline; p), ComputeZ(CF; p)
    L_B, L_CF = ComputeL(Baseline; p), ComputeL(CF; p)

    κ = (
        Lᵈ = log.(CF.Lᵈ[us, 2:K + 1] ./ Baseline.Lᵈ[us, 2:K + 1]) ./ σ_p,
        Lᶠ = log.(CF.Lᶠ[us, 2:K + 1] ./ Baseline.Lᶠ[us, 2:K + 1]) ./ σ_p,
        Wᵈ = log.(CF.Wᵈ[us, 2:K + 1] ./ Baseline.Wᵈ[us, 2:K + 1]) ./ σ_p,
        Wᶠ = log.(CF.Wᶠ[us, 2:K + 1] ./ Baseline.Wᶠ[us, 2:K + 1]) ./ σ_p,
        Z  = log.(Z_CF[us, 2:K + 1]  ./ Z_B[us, 2:K + 1])  ./ σ_p,
        L  = log.(L_CF[us, 2:K + 1]  ./ L_B[us, 2:K + 1])  ./ σ_p,
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
    for σ_p in probes, outcome in (:Lᵈ, :Lᶠ, :Wᵈ, :Wᶠ, :Z, :L)

        d   = maximum(abs.(kernels[σ_p][outcome] .- κ_ref[outcome]))
        rel = d / maximum(abs.(κ_ref[outcome]))

        antisym = haskey(kernels, -σ_p) ?
            maximum(abs.(kernels[σ_p][outcome] .+ kernels[-σ_p][outcome])) : missing

        push!(rows, (; σ_p, outcome, sup_dist = d, rel_dist = rel, antisym_dist = antisym))

    end

    return kernels, DataFrame(rows)

end

#================================================================
                        SIMULATION (BKM CONVOLUTION)
================================================================#
"""
Draw a length-Tsim vector of i.i.d. standard-normal innovations eₜ -- one shared national
draw per period, since mₜ is a single US-specific mobility cost, not state-specific; every
state sees the same shock and differs only through its kernel's exposure. Matches the eₜ in
the paper's law of motion mₜ₊₁=κ+ψ(mₜ-κ)+σeₜ. Draw once per grid search and reuse the same e
across every (σ,ψ) point, so comparisons across the grid aren't contaminated by different
random draws.
"""
draw_innovations(Tsim::Integer; seed::Integer = 1) = (Random.seed!(seed); randn(Tsim))

"""
S×(burnin+Tsim) matrix of i.i.d. N(0,1) innovations for the S-replication indirect-inference
loop: row s is replication s's full innovation history (oldest first, burn-in included).
Drawn ONCE and reused verbatim at every (ψ,σ) grid point queried -- re-drawing per grid
point would make the objective surface a step function of simulation noise rather than of
(ψ,σ), which is the entire point of "keeping the drawn innovations fixed" in indirect
inference. Draw at the full intended S and take `@view E[1:S_proto, :]` to prototype on a
subset, rather than drawing a separate smaller matrix -- that makes the prototype a literal
subset of the production run (same realizations), so S-convergence can be checked directly
by comparing subset vs. full results instead of by re-running with fresh randomness.
"""
function draw_innovations(S::Integer, Tsim::Integer; burnin::Integer, seed::Integer = 1)
    Random.seed!(seed)
    return randn(S, burnin + Tsim)
end

"""
Convolve a single outcome's kernel (L×K, per-unit-σ response from bkm_kernel) against a drawn
innovation sequence e to build a simulated L×Tsim log-deviation panel, for structural
innovation SD σ. Linear MA(∞) representation truncated at the kernel's horizon K:
y[l,τ] = Σₛ κ[l,s+1]·σ·e[τ-s] for s=0,…,min(τ-1,K-1) -- κ(s) already encodes the full
ψ-driven decay of a single shock's effect on mₜ (see mit_shock), so no further AR(1)
rescaling happens here. This convolution is what replaces re-solving the nonlinear model
per draw: κ depends on ψ and must be rebuilt (a fresh bkm_kernel call) if ψ changes, but for
a fixed κ, rescaling by σ and convolving against e is cheap.
"""
function simulate_panel(κ::AbstractMatrix, e::AbstractVector; σ::Real)

    L, K  = size(κ)
    Tsim  = length(e)
    y     = zeros(L, Tsim)

    for τ in 1:Tsim, s in 0:min(τ - 1, K - 1)
        y[:, τ] .+= κ[:, s + 1] .* σ .* e[τ - s]
    end

    return y

end

"""
Simulate a chosen subset of outcomes in a bkm_kernel NamedTuple against the same shared
innovation draw e -- one national shock hits all outcomes/states simultaneously each period;
only the exposure (κ) differs by outcome and state. `outcomes` defaults to every field in
κ, but the S-replication loop only ever needs the primitives (Lᵈ, Lᶠ, Wᵈ, Wᶠ) that
construct_aggregates consumes -- convolving the diagnostic Z/L kernels too would be a third
of the work thrown away, multiplied by every replication and every grid point.
"""
simulate_panel(κ::NamedTuple, e::AbstractVector; σ::Real, outcomes = propertynames(κ)) =
    NamedTuple(outcome => simulate_panel(getfield(κ, outcome), e; σ) for outcome in outcomes)

#================================================================
                        AGGREGATE CONSTRUCTION FROM SIMULATED PRIMITIVES
================================================================#
"""
Construct λ, Z, L on the simulated panel by applying TaskAggregates_LN/LaborAggregate to
simulated LEVELS of the primitives (Lᵈ, Lᶠ, Wᵈ, Wᶠ) -- exactly what ComputeZ/ComputeL do to a
solved Soln's fields -- rather than separately convolving log-linearized Z and L kernels. Z
and L are nonlinear functions of the primitives; convolving two independently-linearized Z/L
kernels would compound two different linearization errors on top of each other, instead of
the one linearization actually licensed by BKM (linearizing the underlying wage/labor-stock
dynamics), with the model's own exact nonlinear aggregator applied afterward -- which is also
how λ, Z, L are constructed from real (not simulated) data: from levels, not from a log-linear
combination of moments.

`primitives` is a NamedTuple of simulated log-deviation panels (L×Tsim each, as returned by
simulate_panel(κ, e; σ) for a κ containing at least Lᵈ, Lᶠ, Wᵈ, Wᶠ). `Baseline` is the flat,
re-anchored steady state the kernel was built around: since it is flat, the level each log
deviation is taken relative to is just Baseline's first column, for every simulated period.
"""
function construct_aggregates(primitives::NamedTuple, Baseline::Soln; p::Parameters)

    us = 1:p.N - 1
    Lᵈ_lvl = Baseline.Lᵈ[us, 1] .* exp.(primitives.Lᵈ)
    Lᶠ_lvl = Baseline.Lᶠ[us, 1] .* exp.(primitives.Lᶠ)
    Wᵈ_lvl = Baseline.Wᵈ[us, 1] .* exp.(primitives.Wᵈ)
    Wᶠ_lvl = Baseline.Wᶠ[us, 1] .* exp.(primitives.Wᶠ)

    ta = TaskAggregates_LN.(p.ρ, p.μ_z, p.ξ_ω, p.ξ_z, Wᵈ_lvl ./ Wᶠ_lvl)
    λ, oneMinusλ, Z = getproperty.(ta, :λ), getproperty.(ta, :oneMinusλ), getproperty.(ta, :Z)
    L = LaborAggregate.(λ, p.ρ, Lᶠ_lvl, Lᵈ_lvl, oneMinusλ)

    return (; λ, Z, L, Lᵈ = Lᵈ_lvl, Lᶠ = Lᶠ_lvl, Wᵈ = Wᵈ_lvl, Wᶠ = Wᶠ_lvl)

end

#================================================================
                        SINGLE-NATIVITY BARTIK-CARD INSTRUMENT
================================================================#
"""
State l's own foreign-born share of its own total labor force at the flat, re-anchored
pre-shock steady state -- the model analogue of Functions.do's fixed 1990 settlement share
s_{Foreign,l,1990} = L^F_{l,1990}/N_{l,1990} (Card 2001's original, undecomposed
settlement instrument; the model has one migrant origin, "Rest of World", so there is no
country-of-origin breakdown to sum over, unlike the paper's real multi-origin Bartik).
Time-invariant by construction (Baseline is flat everywhere), matching how the real
instrument freezes shares before the sample starts. Not degenerate across states -- this is
each state's own foreign-born intensity, and its cross-sectional dispersion is the entire
source of identifying variation once the year fixed effect absorbs the common national
shift (see bartik_instrument).
"""
function bartik_shares(Baseline::Soln; p::Parameters)

    us = 1:p.N - 1
    Lᵈ, Lᶠ = Baseline.Lᵈ[us, 1], Baseline.Lᶠ[us, 1]
    return Lᶠ ./ (Lᵈ .+ Lᶠ)

end

"""
National aggregate foreign-born growth on a simulated (LEVELS) panel:
g[t] = (Σₗ Lᶠ[l,t] - Σₗ Lᶠ[l,t-1]) / Σₗ Lᶠ[l,t-1], for t = 2,…,Tsim (g[1] = NaN, no lag).
Mirrors Functions.do's fg_agg_<region>, collapsed to a single region ("Foreign", i.e. no
per-origin-country decomposition).
"""
function aggregate_foreign_growth(Lᶠ::AbstractMatrix)

    Tsim  = size(Lᶠ, 2)
    total = vec(sum(Lᶠ, dims = 1))
    g     = fill(NaN, Tsim)
    g[2:end] = (total[2:end] .- total[1:end - 1]) ./ total[1:end - 1]

    return g

end

"""
Leave-one-out variant: state l's national shift excludes l's own stock from both the
numerator and the lagged denominator (Functions.do's Supply_Agg_<region>_LOO), giving an
L×Tsim matrix -- state-specific under LOO, unlike aggregate_foreign_growth's shared Tsim-
vector. With only one migrant origin, a single large state is a non-trivial share of the
whole national aggregate, so LOO plausibly matters MORE here than in the real multi-origin
instrument (where it's diluted across ~10 origin groups) -- see bartik_instrument's `loo`
keyword.
"""
function aggregate_foreign_growth_loo(Lᶠ::AbstractMatrix)

    L, Tsim = size(Lᶠ)
    total   = vec(sum(Lᶠ, dims = 1))
    g       = fill(NaN, L, Tsim)

    for l in 1:L
        total_loo   = total .- Lᶠ[l, :]
        g[l, 2:end] = (total_loo[2:end] .- total_loo[1:end - 1]) ./ total_loo[1:end - 1]
    end

    return g

end

"""
Bartik[l,t] = sₗ·gₜ (loo=false) or sₗ·g⁻ˡₜ (loo=true) -- the single-origin collapse of
Functions.do's `egen Bartik_1990 = rowtotal(...)` over origin groups. Column 1 is NaN (no
lag for the growth rate). Default is non-LOO, since IRFEstimates.dta -- the real target
moments -- stores only the non-LOO Iv1990 estimates (MakeIRF.do never saves Iv1990_LOO).

Identifying variation survives the year fixed effect because sₗ·gₜ is a genuine
interaction, not additively separable in l and t: the FE removes gₜ's common level but not
each state's differential (sₗ-weighted) exposure to it.
"""
function bartik_instrument(Lᶠ::AbstractMatrix, s::AbstractVector; loo::Bool = false)

    if loo
        return s .* aggregate_foreign_growth_loo(Lᶠ)
    else
        g = aggregate_foreign_growth(Lᶠ)
        return s .* g'
    end

end

#================================================================
                        SIMULATED REGRESSION PANEL
================================================================#
"""
Assemble the regression-ready simulated panel from construct_aggregates' LEVELS output.

  N[l,t]      = Lᵈ[l,t] + Lᶠ[l,t]                          (Supply_Total analogue)
  fg[l,t]     = (Lᶠ[l,t] - Lᶠ[l,t-1]) / N[l,t-1]           (t=1 -> NaN; matches
                PreRegProcessing's `fg = (f-f[_n-1])/emp[_n-1]` EXACTLY -- a flow-over-
                LAGGED-total-employment ratio, NOT a log difference; do not confuse with
                emp below, which uses CONTEMPORANEOUS N)
  Bartik[l,t] = bartik_instrument(Lᶠ, s; loo)               (t=1 -> NaN)
  emp[l,t]    = N[l,t]                                       (CONTEMPORANEOUS, matching
                the regression weight [pw=emp] -- deliberately a different N-timing than
                fg's denominator above)
  Y           = NamedTuple of the requested outcome LEVEL matrices (from `agg`)
"""
function simulated_panel(agg::NamedTuple, s::AbstractVector; loo::Bool = false, outcomes = (:Z, :L))

    Lᵈ, Lᶠ = agg.Lᵈ, agg.Lᶠ
    N = Lᵈ .+ Lᶠ
    L, Tsim = size(Lᶠ)

    fg = fill(NaN, L, Tsim)
    fg[:, 2:end] = (Lᶠ[:, 2:end] .- Lᶠ[:, 1:end - 1]) ./ N[:, 1:end - 1]

    Bartik = bartik_instrument(Lᶠ, s; loo)
    Y = NamedTuple(outcome => getfield(agg, outcome) for outcome in outcomes)

    return (; fg, Bartik, emp = N, Y)

end

#================================================================
                        WEIGHTED 2SLS / LOCAL PROJECTION
================================================================#
"""
Weighted within-year demeaning of every column of M in place, over the sample rows given:
  x̃[i] = x[i] - (Σⱼ:yr[j]=yr[i] w[j]x[j]) / (Σⱼ:yr[j]=yr[i] w[j]).
By Frisch-Waugh-Lovell this is exactly equivalent to including i.year dummies in the
weighted regression -- since every column of M (LHS, endogenous regressor, instruments,
controls) is demeaned identically, the fixed effect is absorbed for every stage at once,
with no dummy columns and no risk of a collinear/degenerate dummy (unlike the discarded
prior-art approach, which pruned degenerate year dummies after the fact).
"""
function absorb_year!(M::AbstractMatrix, yr::AbstractVector, w::AbstractVector)

    for y in unique(yr)
        idx  = findall(==(y), yr)
        wsum = sum(@view w[idx])
        for j in 1:size(M, 2)
            wmean = sum(@view(w[idx]) .* @view(M[idx, j])) / wsum
            @view(M[idx, j]) .-= wmean
        end
    end

    return M

end

"""
Build the regression-ready long-format arrays for one horizon h from wide L×Tsim panel
matrices, replicating EstimateIRF's sample/lag structure in closed form -- no DataFrame
join (the discarded prior-art `add_shift` used a leftjoin per lag, which is both
unnecessary here since the panel is balanced/gap-free by construction, and far too slow at
S×outcomes×horizons repetitions).

The estimation window is t ∈ [Lmax+2, Tsim-h], where Lmax = max(depvarlags, fglags,
ivlags): the smallest t at which every lag AND the t+h lead are simultaneously defined.
This is an EXACT closed-form match to how Stata's L(k). operators (combined with the real
panel's own absent pre-1994 rows) trim `ivreg2`'s sample -- see the validation checklist
in SMM_Grid.jl for the corresponding per-horizon sample-size check.

  y      = Δ_h ln Y[l,t]  = log(Y[l,t+h]/Y[l,t-1])
  fg     = fg[l,t]                                          (endogenous)
  Bartik = Bartik[l,t]                                       (excluded instrument,
                                                              CONTEMPORANEOUS only)
  Xd     = [L(1:depvarlags).D0Y]     where D0Y = log(Y/L.Y)  (log first-difference lags)
  Xf     = [L(1:fglags).fg]                                   (fg's own lags, in levels)
  Xb     = [L(1:ivlags).Bartik]                                (Bartik's own lags, in
                                                              levels -- exogenous CONTROLS,
                                                              not additional instruments,
                                                              matching EstimateIRF's
                                                              `exogenous(L(1/ivlags).Bartik)`)
  wt     = emp[l,t]                                            (contemporaneous weight)
  cl     = state index l                                        (cluster id)
  yr     = t                                                     (fixed-effect id)
"""
function lpiv_sample(Ylvl::AbstractMatrix, panel::NamedTuple; h::Integer,
    depvarlags::Integer, fglags::Integer, ivlags::Integer)

    L, Tsim = size(Ylvl)
    Lmax    = max(depvarlags, fglags, ivlags)
    trange  = (Lmax + 2):(Tsim - h)

    D0Y = fill(NaN, L, Tsim)
    D0Y[:, 2:end] = log.(Ylvl[:, 2:end] ./ Ylvl[:, 1:end - 1])

    n  = L * length(trange)
    y  = Vector{Float64}(undef, n)
    fg = Vector{Float64}(undef, n)
    Bt = Vector{Float64}(undef, n)
    Xd = Matrix{Float64}(undef, n, depvarlags)
    Xf = Matrix{Float64}(undef, n, fglags)
    Xb = Matrix{Float64}(undef, n, ivlags)
    wt = Vector{Float64}(undef, n)
    cl = Vector{Int}(undef, n)
    yr = Vector{Int}(undef, n)

    i = 0
    for t in trange, l in 1:L
        i += 1
        y[i]  = log(Ylvl[l, t + h] / Ylvl[l, t - 1])
        fg[i] = panel.fg[l, t]
        Bt[i] = panel.Bartik[l, t]
        for k in 1:depvarlags; Xd[i, k] = D0Y[l, t - k];        end
        for k in 1:fglags;     Xf[i, k] = panel.fg[l, t - k];   end
        for k in 1:ivlags;     Xb[i, k] = panel.Bartik[l, t - k]; end
        wt[i] = panel.emp[l, t]
        cl[i] = l
        yr[i] = t
    end

    return (; y, fg, Bartik = Bt, Xd, Xf, Xb, wt, cl, yr)

end

"""
First-stage Wald F for the excluded instrument(s) W_excl in a regression of x on
[W_excl, X_exog], weighted by w with an optional cluster-robust sandwich. Reduces to the
squared t-stat on Bartik when there is exactly one excluded instrument (this build's case:
fg is instrumented by contemporaneous Bartik only, with Bartik's own lags entering as
exogenous controls instead -- see lpiv_sample). Kept SEPARATE from iv_2sls's main
clustered-SE computation because it operates on x as the dependent variable, not y.
"""
function first_stage_F(x::AbstractVector, X_exog::AbstractMatrix, W_excl::AbstractMatrix,
    w::AbstractVector; clusters = nothing)

    n = length(x)
    sw = sqrt.(w)
    Wfull = hcat(W_excl, X_exog)
    Wt = sw .* Wfull
    xt = sw .* x

    π̂     = qr(Wt) \ xt
    resid = xt .- Wt * π̂
    k     = size(Wt, 2)
    n_instr = size(W_excl, 2)

    WtW_inv = inv(Wt' * Wt)

    if isnothing(clusters)
        σ̂² = sum(abs2, resid) / (n - k)
        V  = σ̂² .* WtW_inv
    else
        meat = zeros(k, k)
        for g in unique(clusters)
            idx = findall(==(g), clusters)
            Wg, ug = @view(Wt[idx, :]), @view(resid[idx])
            meat .+= (Wg' * ug) * (Wg' * ug)'
        end
        G   = length(unique(clusters))
        dfc = (G / (G - 1)) * ((n - 1) / (n - k))
        V   = dfc .* (WtW_inv * meat * WtW_inv)
    end

    b  = π̂[1:n_instr]
    Vb = V[1:n_instr, 1:n_instr]

    return (b' * (Vb \ b)) / n_instr

end

"""
Weighted 2SLS via Hansen's projection formula, solved by QR rather than normal equations
(the simulated design is far more collinear than real data -- every regressor is a
different linear combination of the same single national shock history, since there is
only one aggregate mₜ shock in this model -- so squaring the condition number via X'X is
worth avoiding; see the pilot conditioning check in SMM_Grid.jl).

    ŷ = X̂β,  X̂ = P_W X = W(W'W)⁻¹W'X,  β = (X̂'X̂)⁻¹X̂'y

Built for exactly ONE endogenous regressor `x` (fg), matching this model -- not a general
multi-endogenous IV routine. `X_exog` are included exogenous controls (entered in both
stages, i.e. also part of the instrument set W). `W_excl` are the excluded instrument(s)
(entered in the first stage only). `w` are levels weights (the √-scaling happens inside).
Residuals for the covariance use the ACTUAL (unprojected) X, per the standard 2SLS
asymptotic variance formula -- only the "bread" uses the projected X̂.

SEs/clustering/F are skipped entirely unless `se=true`: in the S-replication loop only β̂
is a moment, and the sandwich (plus first_stage_F, itself a second regression) is a
non-trivial fraction of the per-regression cost at S×outcomes×horizons repetitions.
"""
function iv_2sls(y::AbstractVector, x::AbstractVector, X_exog::AbstractMatrix,
    W_excl::AbstractMatrix, w::AbstractVector; clusters = nothing, se::Bool = false)

    n = length(y)
    sw = sqrt.(w)

    Xfull = hcat(x, X_exog)
    Wfull = hcat(W_excl, X_exog)

    yt = sw .* y
    Xt = sw .* Xfull
    Wt = sw .* Wfull

    X̂t = Wt * (qr(Wt) \ Xt)
    β  = qr(X̂t) \ yt

    result = (; β = β[1])
    !se && return result

    û = yt .- Xt * β
    k = size(X̂t, 2)
    X̂tX̂t_inv = inv(X̂t' * X̂t)

    if isnothing(clusters)
        σ̂² = sum(abs2, û) / (n - k)
        V  = σ̂² .* X̂tX̂t_inv
    else
        meat = zeros(k, k)
        for g in unique(clusters)
            idx = findall(==(g), clusters)
            Xg, ug = @view(X̂t[idx, :]), @view(û[idx])
            meat .+= (Xg' * ug) * (Xg' * ug)'
        end
        G   = length(unique(clusters))
        dfc = (G / (G - 1)) * ((n - 1) / (n - k))
        V   = dfc .* (X̂tX̂t_inv * meat * X̂tX̂t_inv)
    end

    F = first_stage_F(x, X_exog, W_excl, w; clusters)

    return merge(result, (; stderr = sqrt(V[1, 1]), V, F))

end

"""
One horizon-h LPIV/2SLS regression on a panel (simulated OR real, via the same wide L×Tsim
matrix convention -- see load_real_panel for the real-data adapter), replicating
EstimateIRF exactly:

  LHS    Δ_h ln Y = log(Y[l,t+h]/Y[l,t-1])                     (cumulative long difference)
  Endog  fg[l,t]                     instrumented by Bartik[l,t] (contemporaneous only)
  Exog   L(1:depvarlags).D0Y, L(1:fglags).fg, L(1:ivlags).Bartik
         (D0Y = log(Y/L.Y) -- log first difference; fg and Bartik lags are the raw
         level/ratio series, NOT logged, matching how they're already-differenced flows)
  FE     year, absorbed by weighted within-year demeaning (Frisch-Waugh-Lovell)
  Wts    contemporaneous emp[l,t] (matches [pw=emp]; weighted=false uses equal weights)
"""
function lpiv_beta(Ylvl::AbstractMatrix, panel::NamedTuple; h::Integer,
    depvarlags::Integer = 3, fglags::Integer = depvarlags, ivlags::Integer = depvarlags,
    weighted::Bool = true, se::Bool = false, ols::Bool = false)

    smp = lpiv_sample(Ylvl, panel; h, depvarlags, fglags, ivlags)
    w   = weighted ? smp.wt : ones(length(smp.y))

    M = hcat(smp.y, smp.fg, smp.Bartik, smp.Xd, smp.Xf, smp.Xb)
    absorb_year!(M, smp.yr, w)

    y      = M[:, 1]
    x      = M[:, 2]
    W_excl = ols ? reshape(x, :, 1) : M[:, 3:3]   # x instrumenting itself == OLS
    X_exog = M[:, 4:end]

    return iv_2sls(y, x, X_exog, W_excl, w; clusters = se ? smp.cl : nothing, se)

end

"""
Full LPIV impulse response, h=0,…,H, on one outcome's (levels) panel -- calls lpiv_beta at
each horizon SEPARATELY (true local projections, not a joint system), matching MakeIRF.do's
per-horizon EstimateIRF calls exactly. Returns the length-(H+1) vector of β coefficients.
"""
lpiv_irf(Ylvl::AbstractMatrix, panel::NamedTuple; H::Integer = 9, kwargs...) =
    [lpiv_beta(Ylvl, panel; h, kwargs...).β for h in 0:H]

#================================================================
                        SIMULATED MOMENTS (S-REPLICATION)
================================================================#
"""
S-averaged simulated LPIV coefficients β̂ʰ(σ,ψ) -- the model moments for indirect
inference. `κ` is a pre-built bkm_kernel for the target ψ (the expensive, ψ-dependent
object -- build ONCE per ψ outside this function, reuse across every σ); `E` is a FIXED
S×(burnin+Tsim) innovation matrix (see draw_innovations), identical at every (σ,ψ) grid
point queried -- only σ enters here as a per-call rescaling.

Per replication s: convolve the primitives, drop the first `burnin` columns (see
simulate_panel/moment_grid for why B=K-1 exactly eliminates the truncated-convolution
cold-start bias), build LEVELS via construct_aggregates, assemble the regression panel via
simulated_panel, then run lpiv_irf per requested outcome. β_reps is kept (not discarded)
since it is exactly the simulation-variance estimate a later SMM weighting matrix needs, at
no extra simulation cost.
"""
function simulated_moments(κ::NamedTuple, Baseline::Soln, E::AbstractMatrix;
    p::Parameters, σ::Real, burnin::Integer, s_share::AbstractVector,
    H::Integer = 9, outcomes = (:Z, :L),
    depvarlags::Integer = 3, fglags::Integer = depvarlags, ivlags::Integer = depvarlags,
    loo::Bool = false, weighted::Bool = true)

    S = size(E, 1)
    β_reps = NamedTuple(o => zeros(S, H + 1) for o in outcomes)
    ok     = trues(S)

    Threads.@threads for s in 1:S

        try
            prims_full = simulate_panel(κ, @view(E[s, :]); σ, outcomes = (:Lᵈ, :Lᶠ, :Wᵈ, :Wᶠ))
            prims = NamedTuple(k => v[:, burnin + 1:end] for (k, v) in pairs(prims_full))
            agg   = construct_aggregates(prims, Baseline; p)
            pnl   = simulated_panel(agg, s_share; loo, outcomes)

            for o in outcomes
                β_reps[o][s, :] = lpiv_irf(getfield(agg, o), pnl; H, depvarlags, fglags, ivlags, weighted)
            end
        catch
            ok[s] = false
        end

    end

    n_ok = count(ok)
    β = NamedTuple(o => vec(sum(β_reps[o][ok, :], dims = 1)) ./ n_ok for o in outcomes)

    return (; β, β_reps, ok, n_ok)

end

"""
(ψ,σ) sensitivity grid of simulated LPIV moments β̂ʰ(σ,ψ). Outer loop over ψ builds ONE
bkm_kernel per ψ (~145s per the last session's timing -- the only expensive step:
SolveSteadyState is never re-run here, since ψ enters only mit_shock and the cached
Baseline is valid for every ψ, see bkm_kernel's ψ keyword); the inner σ loop reuses that
kernel (cheap rescaling + reconvolution only, via simulated_moments).

Returns a long DataFrame (ψ, σ, outcome, h, β, β_sd) -- β_sd is the across-replication
standard deviation at that horizon, a simulation-noise diagnostic (how large is S buying
you), not a structural standard error.
"""
function moment_grid(Baseline::Soln, E::AbstractMatrix; p::Parameters,
    ψ_grid, σ_grid, K::Integer = 60, σ_p::Real = 0.01, burnin::Integer = K - 1,
    CF_tol::Real = 1e-6, CF_maxiter::Integer = 20_000, verbose::Bool = false, kwargs...)

    s_share = bartik_shares(Baseline; p)
    rows = NamedTuple[]

    for ψ in ψ_grid

        κ, _ = bkm_kernel(Baseline; p, σ_p, K, ψ, CF_tol, CF_maxiter, verbose)

        for σ in σ_grid

            mom = simulated_moments(κ, Baseline, E; p, σ, burnin, s_share, kwargs...)
            for o in keys(mom.β), h in 0:length(mom.β[o]) - 1
                push!(rows, (; ψ, σ, outcome = o, h,
                    β    = mom.β[o][h + 1],
                    β_sd = std(@view mom.β_reps[o][mom.ok, h + 1])))
            end

        end

    end

    return DataFrame(rows)

end

#================================================================
                        REAL-DATA VALIDATION ADAPTER
================================================================#
"""
Load StateAnalysisRegReady.dta (PreRegProcessing's saved output, matches the real empirical
regression exactly) and reshape it into the same wide L×Tsim (state × year) matrix
convention lpiv_irf/simulated_panel use, so the IDENTICAL Julia estimator can be run on
REAL data and checked against IRFEstimates.dta's stored (ivreg2-produced) β/SE/F -- the
highest-value single correctness check available for the hand-rolled 2SLS, independent of
any question about the simulation itself. Assumes the panel is balanced (confirmed: 1428
rows = 51 states × 28 years, 1994-2021, no gaps), so row-adjacency lags coincide exactly
with Stata's xtset-based L(k). operators.
"""
function load_real_panel(; outcome::Symbol = :Z)

    outcome in (:Z, :L) || error("outcome must be :Z or :L")

    df = DataFrame(load(joinpath(data, "StateAnalysisRegReady.dta")))
    sort!(df, [:state, :year])

    states = sort(unique(df.state))
    years  = sort(unique(df.year))
    L, Tsim = length(states), length(years)
    sidx = Dict(s => i for (i, s) in enumerate(states))
    tidx = Dict(y => i for (i, y) in enumerate(years))

    Ylvl   = fill(NaN, L, Tsim)
    fg     = fill(NaN, L, Tsim)
    Bartik = fill(NaN, L, Tsim)
    emp    = fill(NaN, L, Tsim)
    ycol   = outcome == :Z ? df.Z : df.L

    for row in 1:nrow(df)
        i, j = sidx[df.state[row]], tidx[df.year[row]]
        Ylvl[i, j]   = coalesce(ycol[row], NaN)
        fg[i, j]     = coalesce(df.fg[row], NaN)
        Bartik[i, j] = coalesce(df.Bartik_1990[row], NaN)
        emp[i, j]    = coalesce(df.emp[row], NaN)
    end

    return Ylvl, (; fg, Bartik, emp), years

end

#================================================================
                SYNTHETIC 2SLS VALIDATION (KNOWN-TRUTH MONTE CARLO)
================================================================#
"""
A from-scratch, closed-form panel DGP -- NOT the CDP/BKM structural model -- built to
validate the hand-rolled `iv_2sls`/`lpiv_beta` machinery against a KNOWN true cumulative IRF,
under a genuine endogeneity problem: the migration-flow regressor `fg` is correlated with the
structural equation's error by construction, exactly the identification problem the LPIV's
two-stage (Bartik-instrumented) structure is meant to solve.

Construction: `fg[l,t] = π₁·Bartik[l,t] + v[l,t]`, where `Bartik[l,t] = s[l]·g[t]` is driven
only by the (exogenous, iid-across-t) national shock `g`. The period-t GROWTH of `Y` responds
to `fg` through a decaying distributed lag `θ_k = θ₀·decay^k` (k=0,…,Kmax), plus a disturbance
`u[l,t] = ρ·v[l,t] + √(1-ρ²)·ε[l,t]` that shares correlation `ρ` with `v` -- so `fg[l,t]` is
correlated with the SAME-period disturbance entering `Y`'s growth (`Cov(fg,u) = ρ·σᵥ² ≠ 0` for
ρ≠0), which is exactly what would bias OLS. `Bartik` itself is a function only of `(s,g)`,
independent of `(v,ε,u)`, so it remains a valid instrument throughout. A common national shock
`η[t]` (loading `δ`) is added to `Y`'s growth too, exercising `absorb_year!`'s year-FE removal.

Because `fg` is iid across t (both `g` and `v` are), the coefficient on `fg[l,t]` in a
regression of `Δ_h ln Y[l,t] = ln Y[l,t+h] - ln Y[l,t-1]` on `fg[l,t]` (plus lags, which have
zero true coefficients here -- there is no serial structure to pick up) is exactly the
cumulative sum `β_h = Σ_{k=0}^{min(h,Kmax)} θ_k`: other-period `fg`'s entering the Δ_h window
are uncorrelated with `fg[l,t]` and only add residual noise, not bias, to that coefficient.
Returns `(Ylvl, panel, β_h_true)` in exactly `lpiv_sample`/`lpiv_irf`'s expected shapes, so
the SAME estimator code used on real and BKM-simulated data runs on this panel unmodified.
"""
function synthetic_iv_panel(s::AbstractVector, Tsim::Integer; seed::Integer,
    H::Integer = 9, θ₀::Real = 0.3, decay::Real = 0.7, Kmax::Integer = 6,
    π₁::Real = 10.0, ρ::Real = 0.6, σ_v::Real = 1.0, σ_ε::Real = 1.0,
    σ_η::Real = 0.5, δ::Real = 1.0, σ_g::Real = 1.0)

    Random.seed!(seed)
    L = length(s)

    g      = σ_g .* randn(Tsim)             # national shock driving the instrument, iid over t
    η      = σ_η .* randn(Tsim)             # national common shock in Y's growth (tests year-FE)
    Bartik = s * g'                         # L×Tsim, independent of (v, ε, u) below

    v  = σ_v .* randn(L, Tsim)              # non-instrument-driven component of fg
    ε  = σ_ε .* randn(L, Tsim)
    u  = ρ .* v .+ sqrt(1 - ρ^2) .* ε       # correlated with v ⟹ correlated with fg: the endogeneity
    fg = π₁ .* Bartik .+ v

    θ = θ₀ .* decay .^ (0:Kmax)

    dlnY = zeros(L, Tsim)
    for t in 1:Tsim, k in 0:min(Kmax, t - 1)
        dlnY[:, t] .+= θ[k + 1] .* fg[:, t - k]
    end
    dlnY .+= u .+ δ .* η'

    Ylvl = ones(L, Tsim)
    for t in 2:Tsim
        Ylvl[:, t] = Ylvl[:, t - 1] .* exp.(dlnY[:, t])
    end

    β_h_true = [sum(θ[1:min(h, Kmax) + 1]) for h in 0:H]
    panel    = (; fg, Bartik, emp = ones(L, Tsim))

    return Ylvl, panel, β_h_true

end
