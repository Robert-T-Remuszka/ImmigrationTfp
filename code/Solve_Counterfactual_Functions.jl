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
Update the counterfactual choice probabilities π̃, given the just-updated hat value changes
(Û := CF.U̇/Baseline.U̇, derived internally) and the fixed, already-solved Baseline. Mirrors
UpdateChoiceProbabilities's vectorized structure, with the baseline's dot probability change
(Πᵈ_dot) and the counterfactual's own lagged probability (Π̃ᵈ_lag) both entering multiplicatively,
matching the paper's π̃ⁿ equation exactly.
"""
function UpdateProbabilitiesHat(CF::Soln, Baseline::Soln, M̂::Vector; p::Parameters)

    (; β, νᵈ, νᶠ, N, Πᵈ₋, Πᶠ₋) = p
    (; Πᵈ, Πᶠ, U̇ᵈ, U̇ᶠ, T) = Baseline
    Π̃ᵈ, Π̃ᶠ = CF.Πᵈ, CF.Πᶠ

    Ûᵈ, Ûᶠ = CF.U̇ᵈ ./ U̇ᵈ, CF.U̇ᶠ ./ U̇ᶠ
    Π̃ᵈ_new, Π̃ᶠ_new = copy(Π̃ᵈ), copy(Π̃ᶠ)

    for t in 1:T - 1

        Πᵈ_lag = t == 1 ? Πᵈ₋ : Πᵈ[:, :, t - 1]
        Πᶠ_lag = t == 1 ? Πᶠ₋ : Πᶠ[:, :, t - 1]
        Π̃ᵈ_lag = t == 1 ? Πᵈ₋ : Π̃ᵈ[:, :, t - 1]
        Π̃ᶠ_lag = t == 1 ? Πᶠ₋ : Π̃ᶠ[:, :, t - 1]

        Πᵈ_dot = Πᵈ[:, :, t] ./ Πᵈ_lag
        Πᶠ_dot = Πᶠ[:, :, t] ./ Πᶠ_lag
        C      = cost_matrix(M̂[t], N)

        Numᵈ = Πᵈ_dot .* Π̃ᵈ_lag .* (Ûᵈ[:, t + 1] .^ (β / νᵈ))' .* C .^ (-1 / νᵈ)
        Numᶠ = Πᶠ_dot .* Π̃ᶠ_lag .* (Ûᶠ[:, t + 1] .^ (β / νᶠ))' .* C .^ (-1 / νᶠ)

        Π̃ᵈ_new[:, :, t] = Numᵈ ./ sum(Numᵈ, dims = 2)
        Π̃ᶠ_new[:, :, t] = Numᶠ ./ sum(Numᶠ, dims = 2)

    end

    return Soln(CF.Wᵈ, CF.Wᶠ, CF.Lᵈ, CF.Lᶠ, CF.U̇ᵈ, CF.U̇ᶠ, Π̃ᵈ_new, Π̃ᶠ_new, CF.T)

end

"""
Update the counterfactual value changes, by backward recursion, given the just-updated
counterfactual probabilities and the fixed Baseline. Returns a Soln whose U̇ᵈ, U̇ᶠ hold
the counterfactual's *actual* value changes (Ũ̇ = Û·U̇), not the hats themselves — Û is always
recovered on demand as CF.U̇/Baseline.U̇, matching UpdateProbabilitiesHat's convention.
"""
function UpdateChangesHat(CF::Soln, Baseline::Soln, M̂::Vector; p::Parameters)

    (; β, νᵈ, νᶠ, N, Πᵈ₋, Πᶠ₋) = p
    (; Πᵈ, Πᶠ, Wᵈ, Wᶠ, U̇ᵈ, U̇ᶠ, T) = Baseline
    Π̃ᵈ, Π̃ᶠ, W̃ᵈ, W̃ᶠ = CF.Πᵈ, CF.Πᶠ, CF.Wᵈ, CF.Wᶠ

    Ûᵈ_new, Ûᶠ_new = ones(N, T), ones(N, T)

    for t in T - 1:-1:1

        Πᵈ_lag = t == 1 ? Πᵈ₋ : Πᵈ[:, :, t - 1]
        Πᶠ_lag = t == 1 ? Πᶠ₋ : Πᶠ[:, :, t - 1]
        Π̃ᵈ_lag = t == 1 ? Πᵈ₋ : Π̃ᵈ[:, :, t - 1]
        Π̃ᶠ_lag = t == 1 ? Πᶠ₋ : Π̃ᶠ[:, :, t - 1]

        Πᵈ_dot = Πᵈ[:, :, t] ./ Πᵈ_lag
        Πᶠ_dot = Πᶠ[:, :, t] ./ Πᶠ_lag
        C      = cost_matrix(M̂[t], N)

        innerᵈ = (Πᵈ_dot .* Π̃ᵈ_lag .* C .^ (-1 / νᵈ)) * (Ûᵈ_new[:, t + 1] .^ (β / νᵈ))
        innerᶠ = (Πᶠ_dot .* Π̃ᶠ_lag .* C .^ (-1 / νᶠ)) * (Ûᶠ_new[:, t + 1] .^ (β / νᶠ))

        ẇᵈ  = t == 1 ? ones(N) : Wᵈ[:, t] ./ Wᵈ[:, t - 1]
        ẇ̃ᵈ = t == 1 ? ones(N) : W̃ᵈ[:, t] ./ W̃ᵈ[:, t - 1]
        ẇᶠ  = t == 1 ? ones(N) : Wᶠ[:, t] ./ Wᶠ[:, t - 1]
        ẇ̃ᶠ = t == 1 ? ones(N) : W̃ᶠ[:, t] ./ W̃ᶠ[:, t - 1]

        ŵᵈ = ẇ̃ᵈ ./ ẇᵈ
        ŵᶠ = ẇ̃ᶠ ./ ẇᶠ

        Ûᵈ_new[:, t] = ŵᵈ .* innerᵈ .^ νᵈ
        Ûᶠ_new[:, t] = ŵᶠ .* innerᶠ .^ νᶠ

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
    err, prev_err, iter, α = 1 + CF_tol, Inf, 0, 0.5

    while err > CF_tol && iter < CF_maxiter

        Ûᵈ_prev, Ûᶠ_prev = CF.U̇ᵈ ./ Baseline.U̇ᵈ, CF.U̇ᶠ ./ Baseline.U̇ᶠ

        CF = UpdateProbabilitiesHat(CF, Baseline, M̂; p)
        CF = UpdateLaborSupply(CF; p)
        CF = SolveTempEq(CF; p, verbose = verbose && iter % 50 == 0)
        CF = UpdateChangesHat(CF, Baseline, M̂; p)

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
