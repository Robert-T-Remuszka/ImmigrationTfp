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

# NOTE (2026-09-14): this file used to also define a self-contained `Parameters`
# struct/constructor here, kept deliberately separate from Solve_Baseline_Functions.jl's
# struct while the baseline/counterfactual solvers still ran on the old Pareto
# parameterization. That wiring step has now happened: Solve_Baseline_Functions.jl's
# `Parameters` struct itself was switched over to (ρ, μ_z, ξ_ω, ξ_z), so a second
# `Parameters` type here would collide (same name, different fields) whenever both
# files are included together. This file now only provides the math library
# (task_args, log_task_integrals, TaskAggregates_LN, LaborAggregate, valid_params)
# consumed by both AggSupply_Estimate.jl (with bare scalars) and
# Solve_Baseline_Functions.jl (via its own `Parameters`/`TaskAggregates` wrapper).
