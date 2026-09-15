#================================================================
Log-normal aggregate-supply block. Replaces the old Pareto/DFS parameterization
(z(τ)=τ^γ, ω Pareto with minimum μ) with Assumption A4 of the paper:
ω̃=ln ω ~ N(0,ξ_ω²), z(τ)=exp(μ_z+ξ_z·Φ⁻¹(τ)). Formulas below are transcribed
from notes/ProdFunc Reparamertize v3.md's boxed Z,λ (equivalently I_1,I_2)
closed-form Gaussian-integral results, with μ_ω≡0 per A4. Verified independently
by re-deriving n=ln(L^D/L^F) from the ratio of the two labor-market-clearing
conditions (paper's Appendix A.2): L^D/L^F = w^{-1/(1-ρ)}·I_2/I_1, so
n = -w̃/(1-ρ) + ln I_2 - ln I_1 (see RelSupplyLog below).

I_1 is the foreign-side integral (mirrors the old I_F_μ), I_2 the domestic-side
(mirrors the old I_D_μ). Identity check: d₁+d₂ == b·ξ always (see task_args).
================================================================#

using StatsFuns: normlogcdf, logsumexp

"""
Shared inner quantities of the log-normal task block. μ_ω ≡ 0 throughout (A4).
    b = ρ/(1-ρ);  ξ = √(ξ_ω²+ξ_z²);  w̃ = ln w
    d₁ = (w̃-μ_z+bξ_ω²)/ξ,  d₂ = (-w̃+μ_z+bξ_z²)/ξ
Identity: d₁+d₂ = bξ always (cheap correctness check on any implementation).
"""
function task_args(ρ, μ_z, ξ_ω, ξ_z, w)

    b  = ρ / (1 - ρ)
    ξ  = sqrt(ξ_ω^2 + ξ_z^2)
    w̃  = log(w)
    d₁ = ( w̃ - μ_z + b * ξ_ω^2) / ξ
    d₂ = (-w̃ + μ_z + b * ξ_z^2) / ξ

    return (; b, ξ, w̃, d₁, d₂)

end

"""
log I₁ (foreign-side) and log I₂ (domestic-side) task integrals. Evaluated in
logs since the e^{½b²ξ²} prefactors overflow for realistic parameter values.
    ln I₁ = ½b²ξ_ω²         + logΦ(d₁)
    ln I₂ = bμ_z + ½b²ξ_z²  + logΦ(d₂)
"""
function log_task_integrals(ρ, μ_z, ξ_ω, ξ_z, w)

    (; b, d₁, d₂) = task_args(ρ, μ_z, ξ_ω, ξ_z, w)
    lnI₁ = 0.5 * b^2 * ξ_ω^2           + normlogcdf(d₁)
    lnI₂ = b * μ_z + 0.5 * b^2 * ξ_z^2 + normlogcdf(d₂)

    return (; lnI₁, lnI₂, b)

end

"""
Task-productivity aggregate Z and foreign-born task share λ, log-normal
parameterization. Z = (I₁+I₂)^(1/b), λ = I₁/(I₁+I₂).
DROP-IN REPLACEMENT for the deleted Pareto TaskAggregates_μ(ρ, γ, μ, w): same
return shape (; Z, λ), same scalar-broadcastable call convention.

Also returns `oneMinusλ`, computed as exp(lnI₂-lnS) rather than 1-λ: for
extreme relative wages at high ρ, λ saturates to exactly 1.0 in double
precision, and naively computing 1-λ from that saturated value loses all
precision (catastrophic cancellation), which then makes LaborAggregate's
(1-λ)^(1-ρ) term jump discontinuously instead of varying smoothly. Same class
of fix as the old Pareto code's stable_expm1_ratio, for a different
sub-expression -- pass `oneMinusλ` to LaborAggregate's keyword argument
wherever precision matters (e.g. inside an optimizer), not the naive `1-λ`.

NOTE: this is deliberately `exp(lnI₂-lnS)`, a direct ratio with no
subtraction at all, NOT `-expm1(lnI₁-lnS)`. That alternative looks equally
cancellation-safe but isn't: once lnI₁-lnS exceeds -eps (i.e. λ has already
saturated to 1.0 in double precision, which happens for ρ≳0.958 at this
model's estimated parameters), expm1 of a value indistinguishable from 0
returns exactly 0, so `-expm1(...)` silently returns exactly 0 too --
identical failure mode to the naive `1-λ`, just one step further down.
`exp(lnI₂-lnS)` has no such cliff: lnI₂-lnS is a large negative number out
there, not a value near 0, so exp of it underflows gracefully toward (but
not jumping discontinuously to) 0.
"""
function TaskAggregates_LN(ρ, μ_z, ξ_ω, ξ_z, w)

    (; lnI₁, lnI₂, b) = log_task_integrals(ρ, μ_z, ξ_ω, ξ_z, w)
    lnS = logsumexp((lnI₁, lnI₂))
    Z   = exp(lnS / b)
    λ   = exp(lnI₁ - lnS)
    oneMinusλ = exp(lnI₂ - lnS)

    return (; Z, λ, oneMinusλ)

end

"""
CES aggregate of foreign and domestic labor supplies given the foreign-born
task share λ. 4-argument call signature unchanged from the old ProdFunc.jl.
Callers that already have an accurate `oneMinusλ` (e.g. from TaskAggregates_LN,
see its docstring) should pass it explicitly as the 5th positional argument
rather than rely on the `1-λ` default, which loses precision once λ saturates
near 1. Positional (not keyword) so this broadcasts correctly element-wise
via `LaborAggregate.(λ, ρ, LF, LD, oneMinusλ)` -- Julia does NOT broadcast
keyword arguments, so a keyword-based version would silently pass the whole
vector to every call instead of one element per call.
"""
LaborAggregate(λ, ρ, LF, LD, oneMinusλ = 1 - λ) =
    (λ^(1 - ρ) * LF^ρ + oneMinusλ^(1 - ρ) * LD^ρ)^(1 / ρ)

"""
Model-implied log relative supply n = ln(L^D/L^F) at relative wage w = w^D/w^F.
This is the "Step 2" estimating equation (notes/EstimationProductionBlock.md):
    n̂(w; ρ, μ_z, ξ_ω, ξ_z) = -w̃/(1-ρ) + ln I₂ - ln I₁
With (μ_z,ξ_ω,ξ_z) fixed from the choice-probability probit (EstimateCp.do),
ρ is the only unknown here -- see AggSupply_Estimate.jl.
"""
function RelSupplyLog(ρ, μ_z, ξ_ω, ξ_z, w)

    (; lnI₁, lnI₂) = log_task_integrals(ρ, μ_z, ξ_ω, ξ_z, w)

    return -log(w) / (1 - ρ) + lnI₂ - lnI₁

end

"""
Parameter-admissibility guard. The log-normal distribution has all moments, so
the old Pareto tail-existence condition (b<1, i.e. ρ<0.5) is GONE -- only
ρ<1, ρ≠0, and positive dispersions remain.
"""
valid_params(ρ, ξ_ω, ξ_z) = (ρ < 1) && !isapprox(ρ, 0) && (ξ_ω > 0) && (ξ_z > 0)

#================================================================
                        PARAMETERS STRUCT
================================================================#
"""
New, self-contained parameters struct for the log-normal aggregate-supply and
migration blocks. This is DELIBERATELY SEPARATE from Solve_Baseline_Functions.jl's
existing `Parameters` struct (still the old Pareto parameterization) -- per
2026-09-13 scope decision, this refactor only builds AggSupply_Functions.jl and
AggSupply_Estimate.jl as new, self-contained files; wiring this struct into the
baseline/counterfactual solvers (replacing the old `Parameters`, editing its ~8
call sites in Solve_Baseline_Functions.jl) is a separate, later step, at which
point the old ProdFunc.jl/ProdFunc_Estimate.jl also get retired.

Unlike Solve_Baseline_Functions.jl's `Parameters`, this struct holds only
STRUCTURAL PARAMETERS -- no geography/initial-condition data (N, wᵈ_row,
wᶠ_row, Πᵈ₋, Πᶠ₋, Lᵈ₀, Lᶠ₀, Y₀), since those aren't outputs of any estimating
equation and belong to the (untouched, for now) baseline-solve machinery.

θ, δ are state-indexed vectors (paired with `statefip`), not scalars: θ_l is
backed out state-by-state from the capital FOC (3.9) using El-Shagi-Yamarik's
depreciation series, not estimated -- see CalibrateTheta.do.

ψ, κ, σₘ (the mₜ AR(1) process mₜ₊₁=κ+ψ(mₜ-κ)+σₘeₜ) are PLACEHOLDERS pending
indirect inference against the section-4 LPIV IRFs (see Estimation_Funcs.jl's
load_irf_estimates()). κ=0 since the dynamic hat algebra differences it out
until the counterfactual.
"""
struct Parameters{T1 <: Real}

    β::T1                      # HH discount rate
    r::T1                      # Capital rental rate, 1/β - 1

    ρ::T1                      # CES parameter between foreign/domestic task aggregates
    μ_z::T1                    # Location of ln z(τ) = μ_z + ξ_z·Φ⁻¹(τ)
    ξ_ω::T1                    # SD of ln ω ~ N(0, ξ_ω²)   [μ_ω ≡ 0, Assumption A4]
    ξ_z::T1                    # SD of the comparative-advantage schedule ln z(τ)

    θ::Vector{T1}              # Capital share by state, (r+δ_l)·K/Y, from the capital FOC
    δ::Vector{T1}              # Capital depreciation rate by state (El-Shagi-Yamarik)
    statefip::Vector{String}   # State FIPS codes indexing θ, δ

    νᵈ::T1                     # Gumbel scale - domestic
    νᶠ::T1                     # Gumbel scale - foreign

    ψ::T1                      # PLACEHOLDER: AR(1) persistence of the mobility cost mₜ
    κ::T1                      # PLACEHOLDER: AR(1) long-run mean of mₜ (=0 until counterfactual)
    σₘ::T1                     # PLACEHOLDER: AR(1) innovation SD of mₜ

end

"""
Construct Parameters, defaulting every estimated/calibrated field to the
current saved output of its estimating equation:
    ρ, μ_z, ξ_ω, ξ_z  <- AggSupply_Estimate.jl's AggSupply.jld2 (load_aggsupply_estimate)
    θ, δ, statefip    <- CalibrateTheta.do's ThetaDelta.dta      (load_theta_delta)
    νᵈ, νᶠ            <- EstimateScaleBetaAcs.do's NuBetaEstimatesAcs.dta (load_scale_estimate)
β, r are calibrated (not estimated); ψ, κ, σₘ are placeholders (see struct docstring).

Requires Estimation_Funcs.jl to already be included (for the load_* functions) --
matches the existing convention where functions files assume their dependencies
are loaded by the calling script, not by each other.
"""
function Parameters(;
    β::Real                  = 0.96,
    r::Real                  = 1 / β - 1,
    aggsupply::NamedTuple    = load_aggsupply_estimate(),
    scale::NamedTuple        = load_scale_estimate(),
    θδ::DataFrame            = load_theta_delta(),
    ρ::Real                  = aggsupply.ρ,
    μ_z::Real                = aggsupply.μ_z,
    ξ_ω::Real                = aggsupply.ξ_ω,
    ξ_z::Real                = aggsupply.ξ_z,
    νᵈ::Real                 = scale.νᵈ,
    νᶠ::Real                 = scale.νᶠ,
    ψ::Real                  = 0.50,
    κ::Real                  = 0.00,
    σₘ::Real                 = 1.00,
    )

    # Stata's plain `gen` (used in EstimateCp.do/EstimateScaleBetaAcs.do) stores
    # Float32, not Float64 -- promote everything to a single common type here
    # rather than relying on Julia to unify mixed Float32/Float64 inputs under
    # one type parameter, which it won't do automatically.
    T = Float64
    statefip = String.(θδ.statefip)
    θ = Vector{T}(θδ.theta_l)
    δ = Vector{T}(θδ.delta_l)

    return Parameters(T(β), T(r), T(ρ), T(μ_z), T(ξ_ω), T(ξ_z), θ, δ, statefip,
                       T(νᵈ), T(νᶠ), T(ψ), T(κ), T(σₘ))

end

"""
Wrapper mirroring the old Solve_Baseline_Functions.jl call convention
(`TaskAggregates(w; p::Parameters)`), so wiring this struct into the solvers
later is a one-line change at each call site.
"""
TaskAggregates(w; p::Parameters) = TaskAggregates_LN(p.ρ, p.μ_z, p.ξ_ω, p.ξ_z, w)
