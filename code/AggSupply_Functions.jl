#================================================================
Log-normal aggregate-supply block.

I_1 is the foreign-side integral, I_2 the domestic-side
================================================================#

using StatsFuns: normlogcdf, logsumexp

"""
Shared inner quantities of the log-normal task block. μ_ω ≡ 0 throughout.
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
parameterization. Z = (I₁+I₂)^(1/b), λ = I₁/(I₁+I₂) = I_1 / S where S = I₁ + I₂.

Also returns `oneMinusλ`, computed as exp(lnI₂-lnS) rather than 1-λ: for
extreme relative wages at high ρ, λ saturates to exactly 1.0 in double
precision, and naively computing 1-λ from that saturated value loses all
precision (catastrophic cancellation), which then makes LaborAggregate's
(1-λ)^(1-ρ) term jump discontinuously instead of varying smoothly.
Pass `oneMinusλ` to LaborAggregate's keyword argument wherever precision 
matters (e.g. inside an optimizer), not the naive `1-λ`.
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
task share λ.
Callers that already have an accurate `oneMinusλ` (e.g. from TaskAggregates_LN,
see its docstring) should pass it explicitly as the 5th positional argument
rather than rely on the `1-λ` default. 
Positional (not keyword) so this broadcasts correctly element-wise
via `LaborAggregate.(λ, ρ, LF, LD, oneMinusλ)`.
"""
LaborAggregate(λ, ρ, LF, LD, oneMinusλ = 1 - λ) =
    (λ^(1 - ρ) * LF^ρ + oneMinusλ^(1 - ρ) * LD^ρ)^(1 / ρ)

"""
Parameter-admissibility guard. The log-normal distribution has all moments, so
the old Pareto tail-existence condition (b<1, i.e. ρ<0.5) is GONE -- only
ρ<1, ρ≠0, and positive dispersions remain.
"""
valid_params(ρ, ξ_ω, ξ_z) = (ρ < 1) && !isapprox(ρ, 0) && (ξ_ω > 0) && (ξ_z > 0)