#================================================================
                STATE-LEVEL UNIT-SHOCK IMPULSE RESPONSE
================================================================#
"""
Solve one counterfactual to a σ-sized MIT shock (mit_shock) and return, for every US
location, the log-deviation impulse response of wᴰ, wᶠ, Z, L (task-aggregate labor), and Lᶠ
over horizons k=0,…,K-1 (k=0 is the period the shock first hits), normalized by σ so the
result is a per-unit ("numerical derivative") response, per Boppart, Krusell, and Mitman
(2018). wᴰ, wᶠ, Z, L are the paper's four target IRF variables (Figure 1); Lᶠ is retained
separately since it's what drives the migration-inflow impulse (fg) in the LPIV.
"""
function compute_state_irf(Baseline::Soln; p::Parameters, σ::Real = 1.0, K::Integer = 60,
    CF_maxiter::Integer = 2000, verbose::Bool = false)

    M̂ = [mit_shock(σ = σ, ψ = p.ψ)(t) for t in 1:(Baseline.T - 1)]
    CF = SolveCounterfactual(M̂; Baseline, p, CF_maxiter, verbose)

    us      = 1:(p.N - 1)
    horizon = 2:(K + 1)   # Julia indices; horizon[1] (index 2) is k=0, the period the shock first hits

    logŵᴰ = log.(CF.Wᵈ[us, horizon] ./ Baseline.Wᵈ[us, horizon]) ./ σ
    logŵᶠ = log.(CF.Wᶠ[us, horizon] ./ Baseline.Wᶠ[us, horizon]) ./ σ
    logL̂ᶠ = log.(CF.Lᶠ[us, horizon] ./ Baseline.Lᶠ[us, horizon]) ./ σ

    Z_Baseline, Z_CF = ComputeZ(Baseline; p), ComputeZ(CF; p)
    logẐ = log.(Z_CF[us, horizon] ./ Z_Baseline[us, horizon]) ./ σ

    L_Baseline, L_CF = ComputeL(Baseline; p), ComputeL(CF; p)
    logL̂ = log.(L_CF[us, horizon] ./ L_Baseline[us, horizon]) ./ σ

    return (wᴰ = logŵᴰ, wᶠ = logŵᶠ, Z = logẐ, L = logL̂, Lᶠ = logL̂ᶠ)

end

#================================================================
                BKM PANEL SIMULATION (IRF CONVOLUTION)
================================================================#
"""
Simulate a state-level panel of log-deviations from baseline by convolving each state's own
unit-shock IRF (from compute_state_irf) with one shared, drawn sequence of standardized
innovations eₜ. Per Boppart, Krusell, and Mitman (2018), the model is treated as linear in
the shock, so the response at simulated period T is
    xₗ(T) = σ Σ_{k=0}^{K-1} IRFₗ(k)·e_{T-k},
i.e. each state's simulated path is its own IRF convolved with the same national innovation
sequence, since the shock here is to the single US-specific mobility cost (not state-specific).
Draws K-1 extra pre-sample innovations so the first simulated periods aren't artificially
truncated by a cold start.
"""
function simulate_panel(IRF::NamedTuple, σ::Real, Tsim::Integer; seed::Union{Integer, Nothing} = nothing)

    isnothing(seed) || Random.seed!(seed)

    N, K = size(IRF.wᴰ)
    e = randn(Tsim + K - 1)   # e[j] corresponds to simulated time t = j - (K - 1)

    wᴰ_sim, wᶠ_sim, Z_sim, L_sim, Lᶠ_sim = (zeros(N, Tsim) for _ in 1:5)

    for l in 1:N, T in 1:Tsim, k in 0:(K - 1)
        e_t = e[T - k + K - 1]
        wᴰ_sim[l, T] += σ * IRF.wᴰ[l, k + 1] * e_t
        wᶠ_sim[l, T] += σ * IRF.wᶠ[l, k + 1] * e_t
        Z_sim[l, T]  += σ * IRF.Z[l, k + 1]  * e_t
        L_sim[l, T]  += σ * IRF.L[l, k + 1]  * e_t
        Lᶠ_sim[l, T] += σ * IRF.Lᶠ[l, k + 1] * e_t
    end

    return (wᴰ = wᴰ_sim, wᶠ = wᶠ_sim, Z = Z_sim, L = L_sim, Lᶠ = Lᶠ_sim), e

end