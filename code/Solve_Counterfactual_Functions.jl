#================================================================
                        MOBILITY-COST SHOCK PATH
================================================================#
"""
Proportional mobility-cost path for a one-time MIT/BKM shock: households learn at t=1 of a
σ-sized innovation to the US-specific mobility cost mₜ (e.g. one standard deviation, σ from
mₜ₊₁=κ+ψ(mₜ-κ)+σeₜ), which then decays per the AR(1) law of motion thereafter (eₜ=0 for
t>1). Since the baseline uses no_shock() (Ṁₜ≡1 always), this proportional-change path is
directly usable as M̂ₜ too (M̂ₜ := Ṁ̃ₜ/Ṁₜ = Ṁ̃ₜ/1 = Ṁ̃ₜ).
"""
mit_shock(; σ::Real, ψ::Real) = t -> t == 1 ? exp(σ) : exp(-(1 - ψ) * ψ^(t - 2) * σ)

#================================================================
        COUNTERFACTUAL DYNAMIC BLOCK
================================================================#
"""
Quantities that depend only on the (fixed) Baseline, precomputed once so the counterfactual's
fixed-point loop doesn't rebuild them on every iteration:

  * `Πᵈ_dot`, `Πᶠ_dot` — the baseline's own probability changes Πₜ/Πₜ₋₁, an N×N division per
    period that was previously redone on every one of the hundreds of counterfactual iterations
    despite being identical each time.
  * `ẇᵈ`, `ẇᶠ` — the baseline's gross wage changes Wₜ/Wₜ₋₁ (column 1 left at 1, matching the
    t = 1 convention in UpdateChangesHat).

Both update functions below accept this as a keyword defaulting to a fresh build, so calling
them standalone still works — it is just slower, which is why SolveCounterfactual builds it once
and passes it in.
"""
function counterfactual_cache(Baseline::Soln; p::Parameters)

    (; N, Πᵈ₋, Πᶠ₋) = p
    (; Πᵈ, Πᶠ, Wᵈ, Wᶠ, T) = Baseline

    Πᵈ_dot, Πᶠ_dot = similar(Πᵈ), similar(Πᶠ)
    for t in 1:T - 1
        Πᵈ_lag = t == 1 ? Πᵈ₋ : @view Πᵈ[:, :, t - 1]
        Πᶠ_lag = t == 1 ? Πᶠ₋ : @view Πᶠ[:, :, t - 1]
        @views Πᵈ_dot[:, :, t] .= Πᵈ[:, :, t] ./ Πᵈ_lag
        @views Πᶠ_dot[:, :, t] .= Πᶠ[:, :, t] ./ Πᶠ_lag
    end

    ẇᵈ, ẇᶠ = ones(N, T), ones(N, T)
    for t in 2:T
        @views ẇᵈ[:, t] .= Wᵈ[:, t] ./ Wᵈ[:, t - 1]
        @views ẇᶠ[:, t] .= Wᶠ[:, t] ./ Wᶠ[:, t - 1]
    end

    return (; Πᵈ_dot, Πᶠ_dot, ẇᵈ, ẇᶠ)

end

"""
Update the counterfactual choice probabilities π̃, given the just-updated hat value changes
(Û := CF.U̇/Baseline.U̇, derived internally) and the fixed, already-solved Baseline. Mirrors
UpdateChoiceProbabilities's structure, with the baseline's dot probability change (Πᵈ_dot,
from the cache) and the counterfactual's own lagged probability (Π̃ᵈ_lag) both entering
multiplicatively, matching the paper's π̃ⁿ equation exactly.
"""
function UpdateProbabilitiesHat(CF::Soln, Baseline::Soln, M̂::Vector; p::Parameters,
    cache = counterfactual_cache(Baseline; p))

    (; β, νᵈ, νᶠ, N, Πᵈ₋, Πᶠ₋) = p
    (; U̇ᵈ, U̇ᶠ, T) = Baseline
    (; Πᵈ_dot, Πᶠ_dot) = cache
    Π̃ᵈ, Π̃ᶠ = CF.Πᵈ, CF.Πᶠ

    Π̃ᵈ_new, Π̃ᶠ_new = copy(Π̃ᵈ), copy(Π̃ᶠ)

    Numᵈ, Numᶠ = zeros(N, N), zeros(N, N)
    fᵈ, fᶠ     = zeros(1, N), zeros(1, N)
    sᵈ, sᶠ     = zeros(N, 1), zeros(N, 1)

    for t in 1:T - 1

        Π̃ᵈ_lag = t == 1 ? Πᵈ₋ : @view Π̃ᵈ[:, :, t - 1]
        Π̃ᶠ_lag = t == 1 ? Πᶠ₋ : @view Π̃ᶠ[:, :, t - 1]

        # Û = CF.U̇ ./ Baseline.U̇, formed one column at a time rather than as a full N×T temporary
        @views fᵈ .= (CF.U̇ᵈ[:, t + 1] ./ U̇ᵈ[:, t + 1])' .^ (β / νᵈ)
        @views fᶠ .= (CF.U̇ᶠ[:, t + 1] ./ U̇ᶠ[:, t + 1])' .^ (β / νᶠ)

        @views Numᵈ .= Πᵈ_dot[:, :, t] .* Π̃ᵈ_lag .* fᵈ
        @views Numᶠ .= Πᶠ_dot[:, :, t] .* Π̃ᶠ_lag .* fᶠ
        apply_cost_change!(Numᵈ, M̂[t], νᵈ, N)
        apply_cost_change!(Numᶠ, M̂[t], νᶠ, N)

        sum!(sᵈ, Numᵈ)
        sum!(sᶠ, Numᶠ)

        @views Π̃ᵈ_new[:, :, t] .= Numᵈ ./ sᵈ
        @views Π̃ᶠ_new[:, :, t] .= Numᶠ ./ sᶠ

    end

    return Soln(CF.Wᵈ, CF.Wᶠ, CF.Lᵈ, CF.Lᶠ, CF.U̇ᵈ, CF.U̇ᶠ, Π̃ᵈ_new, Π̃ᶠ_new, CF.T)

end

"""
Update the counterfactual value changes, by backward recursion, given the just-updated
counterfactual probabilities and the fixed Baseline. Returns a Soln whose U̇ᵈ, U̇ᶠ hold
the counterfactual's *actual* value changes (Ũ̇ = Û·U̇), not the hats themselves — Û is always
recovered on demand as CF.U̇/Baseline.U̇, matching UpdateProbabilitiesHat's convention.
"""
function UpdateChangesHat(CF::Soln, Baseline::Soln, M̂::Vector; p::Parameters,
    cache = counterfactual_cache(Baseline; p))

    (; β, νᵈ, νᶠ, N, Πᵈ₋, Πᶠ₋) = p
    (; U̇ᵈ, U̇ᶠ, T) = Baseline
    (; Πᵈ_dot, Πᶠ_dot, ẇᵈ, ẇᶠ) = cache
    Π̃ᵈ, Π̃ᶠ, W̃ᵈ, W̃ᶠ = CF.Πᵈ, CF.Πᶠ, CF.Wᵈ, CF.Wᶠ

    Ûᵈ_new, Ûᶠ_new = ones(N, T), ones(N, T)

    Aᵈ, Aᶠ         = zeros(N, N), zeros(N, N)
    fᵈ, fᶠ         = zeros(N), zeros(N)
    innerᵈ, innerᶠ = zeros(N), zeros(N)

    for t in T - 1:-1:1

        Π̃ᵈ_lag = t == 1 ? Πᵈ₋ : @view Π̃ᵈ[:, :, t - 1]
        Π̃ᶠ_lag = t == 1 ? Πᶠ₋ : @view Π̃ᶠ[:, :, t - 1]

        @views fᵈ .= Ûᵈ_new[:, t + 1] .^ (β / νᵈ)
        @views fᶠ .= Ûᶠ_new[:, t + 1] .^ (β / νᶠ)

        @views Aᵈ .= Πᵈ_dot[:, :, t] .* Π̃ᵈ_lag
        @views Aᶠ .= Πᶠ_dot[:, :, t] .* Π̃ᶠ_lag

        # the inner sum with all costs at 1, then the one row whose costs actually change
        mul!(innerᵈ, Aᵈ, fᵈ)
        mul!(innerᶠ, Aᶠ, fᶠ)
        innerᵈ[N] = inner_row_N(Aᵈ, fᵈ, M̂[t], νᵈ, N)
        innerᶠ[N] = inner_row_N(Aᶠ, fᶠ, M̂[t], νᶠ, N)

        @inbounds for i in 1:N
            # ŵ = (counterfactual gross wage change) / (baseline's, precomputed); both are 1 at t = 1
            ŵᵈ = t == 1 ? one(eltype(W̃ᵈ)) : (W̃ᵈ[i, t] / W̃ᵈ[i, t - 1]) / ẇᵈ[i, t]
            ŵᶠ = t == 1 ? one(eltype(W̃ᶠ)) : (W̃ᶠ[i, t] / W̃ᶠ[i, t - 1]) / ẇᶠ[i, t]
            Ûᵈ_new[i, t] = ŵᵈ * innerᵈ[i]^νᵈ
            Ûᶠ_new[i, t] = ŵᶠ * innerᶠ[i]^νᶠ
        end

    end

    return Soln(CF.Wᵈ, CF.Wᶠ, CF.Lᵈ, CF.Lᶠ, Ûᵈ_new .* U̇ᵈ, Ûᶠ_new .* U̇ᶠ, CF.Πᵈ, CF.Πᶠ, CF.T)

end

#================================================================
                        OUTER COUNTERFACTUAL SOLVER
================================================================#
"""
Solve for the counterfactual sequential equilibrium relative to an already-solved Baseline,
given a path for M̂ (the counterfactual-relative-to-baseline mobility-cost hat — see
mit_shock). The static block needs no new code: UpdateLaborSupply and SolveTempEq are reused
directly, unchanged, since CF shares Baseline's exact t=0 anchor (see the "Why the
Counterfactual Static Block Needs No New Code" note). Only the migration/dynamic block
(UpdateProbabilitiesHat/UpdateChangesHat) is new. Mirrors SolveBaseline's damped
fixed-point structure, damping on Û directly.
"""
function SolveCounterfactual(M̂::Vector; Baseline::Soln, p::Parameters,
    CF_tol::Real = 1e-4, CF_maxiter::Integer = 10_000, verbose::Bool = true)

    CF = deepcopy(Baseline)
    cache = counterfactual_cache(Baseline; p)
    err, prev_err, iter, α = 1 + CF_tol, Inf, 0, 0.5

    while err > CF_tol && iter < CF_maxiter

        Ûᵈ_prev, Ûᶠ_prev = CF.U̇ᵈ ./ Baseline.U̇ᵈ, CF.U̇ᶠ ./ Baseline.U̇ᶠ

        CF = UpdateProbabilitiesHat(CF, Baseline, M̂; p, cache)
        CF = UpdateLaborSupply(CF; p)
        CF = SolveTempEq(CF; p, verbose = verbose && iter % 50 == 0)
        CF = UpdateChangesHat(CF, Baseline, M̂; p, cache)

        Ûᵈ_new, Ûᶠ_new = CF.U̇ᵈ ./ Baseline.U̇ᵈ, CF.U̇ᶠ ./ Baseline.U̇ᶠ
        err   = max(maximum(abs.(Ûᵈ_new .- Ûᵈ_prev)), maximum(abs.(Ûᶠ_new .- Ûᶠ_prev)))
        α_new = err > prev_err ? max(α * 0.5, 0.01) : min(α * 1.01, 0.9)

        Ûᵈ_damped = α_new .* Ûᵈ_new .+ (1 - α_new) .* Ûᵈ_prev
        Ûᶠ_damped = α_new .* Ûᶠ_new .+ (1 - α_new) .* Ûᶠ_prev
        CF = Soln(CF.Wᵈ, CF.Wᶠ, CF.Lᵈ, CF.Lᶠ, Ûᵈ_damped .* Baseline.U̇ᵈ, Ûᶠ_damped .* Baseline.U̇ᶠ, CF.Πᵈ, CF.Πᶠ, CF.T)

        prev_err = err
        α        = α_new
        iter    += 1

        verbose && iter % 50 == 0 && println("****** Counterfactual iteration $iter: err = $err (α = $α) ******")

    end

    err <= CF_tol || error("Counterfactual failed to converge within $CF_maxiter iterations (err = $err) — Baseline's horizon T may need to be longer")

    return CF

end