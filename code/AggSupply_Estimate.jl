# %% Setup
using DataFrames, StatFiles, TidierData, Statistics, ForwardDiff, Optim, JLD2, CSV

include("Globals.jl")
include("Estimation_Funcs.jl")
include("AggSupply_Functions.jl")

#=================================================================
Estimates rho (the CES parameter between the domestic and foreign task
aggregates) from the model's factor-share condition,
    s^F_lt = w^F_lt*L^F_lt / (w^F_lt*L^F_lt + w^D_lt*L^D_lt) = lambda^(1-rho)*(L^F_lt/L_lt)^rho,
(mu_z, xi_omega, xi_z) are fixed/known in advance from EstimateCp.do, so rho
is the only parameter estimated here.

Estimation uses StateWageSupplyPanel.dta (ACS 2001-2024, FTFY, unresidualized
hourly wages). Both sides of the factor-share condition depend on wages only
through the relative wage w=W^D/W^F (RHS, via I^F/I^D) and total labor income
(LHS, w*L -- invariant to how it's split into price x quantity), so the
choice between hourly and annual-equivalent wage units doesn't affect rho.
=================================================================#

# %% 1. Fixed parameters from the prior estimating equation
cp = load_cp_estimate()
println("Fixed from EstimateCp.do: μ_z=", cp.μ_z, " ξ_ω=", cp.ξ_ω, " ξ_z=", cp.ξ_z)

# %% 2. Load and clean the panel -- only wages and labor supplies, nothing else
StateAnalysis = @chain DataFrame(load(joinpath(data, "StateWageSupplyPanel.dta"))) begin
    @mutate(
        Supply_Foreign  = Float64.(Supply_Foreign),
        Supply_Domestic = Float64.(Supply_Domestic),
        Wage_Domestic   = Float64.(Wage_Domestic_Unresid),
        Wage_Foreign    = Float64.(Wage_Foreign_Unresid)
    )
    @rename(statefip = STATEFIP, year = Year)
end

LF = StateAnalysis.Supply_Foreign
LD = StateAnalysis.Supply_Domestic
wF = StateAnalysis.Wage_Foreign
wD = StateAnalysis.Wage_Domestic
w  = wD ./ wF

sF_data = (wF .* LF) ./ (wF .* LF .+ wD .* LD)

# Print up a little summary of the data for the user
println("N obs = ", length(w), "  N states = ", length(unique(StateAnalysis.statefip)))
println("relative wage w: min=", minimum(w), " max=", maximum(w))
println("observed s^F: min=", minimum(sF_data), " max=", maximum(sF_data), " mean=", mean(sF_data))

# %% 3. Model-implied s^F at a candidate ρ, and the NLS objective. No intercept --
#         this is a level (share) equation, not a log-output equation, so there's
#         nothing to profile out.
"""Model-implied factor share ŝ^F = λ^(1-ρ)·(L^F/L)^ρ at every panel row, for a
candidate ρ, with (μ_z,ξ_ω,ξ_z) fixed. Uses oneMinusλ (not 1-λ) throughout --
see AggSupply_Functions.jl's TaskAggregates_LN docstring for why."""
function sF_hat(ρ)
    ta = TaskAggregates_LN.(ρ, cp.μ_z, cp.ξ_ω, cp.ξ_z, w)
    λ  = getproperty.(ta, :λ)
    oneMinusλ = getproperty.(ta, :oneMinusλ)
    L  = LaborAggregate.(λ, ρ, LF, LD, oneMinusλ)
    return λ.^(1 - ρ) .* (LF ./ L).^ρ
end

ssr(ρ) = sum(abs2, sF_data .- sF_hat(ρ))

# %% 4. Grid scan first (sanity check for an interior optimum before trusting Brent)
println()
println("=== SSR grid ===")
for ρ in 0.05:0.05:0.95
    println("ρ=", round(ρ, digits = 2), "  SSR=", round(ssr(ρ), digits = 4))
end

# %% 5. NLS for ρ (Brent, 1-D). ρ ∈ (0,1) the model's own maintained domain
fit = Optim.optimize(ssr, 0.01, 0.99, Optim.Brent())
ρ̂   = Optim.minimizer(fit)

println()
println("Converged: ", Optim.converged(fit))
println("ρ̂ = ", ρ̂, "   SSR = ", Optim.minimum(fit))

# %% 6. Cluster-robust SE for ρ̂ (by state), Gauss-Newton/M-estimation sandwich.
# Score ψ_i = e_i(ρ̂)·h_i(ρ̂), h_i = ∂ŝ^F_i/∂ρ. Hessian approximated by -Σh_i²
# (standard Gauss-Newton, drops residual curvature).
let
    e = sF_data .- sF_hat(ρ̂)
    h = ForwardDiff.derivative(sF_hat, ρ̂)
    ψ = e .* h

    clustered = combine(groupby(DataFrame(statefip = StateAnalysis.statefip, ψ = ψ),
                                 :statefip), :ψ => sum => :ψsum)
    G = nrow(clustered)
    meat  = sum(abs2, clustered.ψsum) * G / (G - 1)
    bread = sum(abs2, h)

    global se_ρ = sqrt(meat) / bread
    println("Clustered SE (", G, " states): se(ρ̂) = ", se_ρ)
end

# %% 7. Diagnostics
println()
println("=== Diagnostics ===")

ta = TaskAggregates_LN.(ρ̂, cp.μ_z, cp.ξ_ω, cp.ξ_z, w)
Z  = getproperty.(ta, :Z)
λ  = getproperty.(ta, :λ)
oneMinusλ = getproperty.(ta, :oneMinusλ)
L  = LaborAggregate.(λ, ρ̂, LF, LD, oneMinusλ)
sF_fit = sF_hat(ρ̂)

println("Z range: ", extrema(Z), "   Z CoV: ", std(Z) / mean(Z))
println("λ range: ", extrema(λ), "   (check: not pinned at exactly 0 or 1)")
println("corr(λ̂, raw foreign share LF/(LF+LD)) = ", cor(λ, LF ./ (LF .+ LD)))
println("corr(ŝ^F, s^F data) = ", cor(sF_fit, sF_data))
println("R² (1 - SSR/TSS) = ", 1 - Optim.minimum(fit) / sum(abs2, sF_data .- mean(sF_data)))

let
    checks = [begin
        (; b, ξ, d₁, d₂) = task_args(ρ̂, cp.μ_z, cp.ξ_ω, cp.ξ_z, wi)
        abs(d₁ + d₂ - b * ξ)
    end for wi in w]
    println("max|d₁+d₂-bξ| = ", maximum(checks), " (should be ~0)")
end

# %% 8. Save to jld2 for future use in model solver
p_star = (ρ = ρ̂, se_ρ = se_ρ,
          μ_z = cp.μ_z, ξ_ω = cp.ξ_ω, ξ_z = cp.ξ_z,
          se_μ_z = cp.se_μ_z, se_ξ_ω = cp.se_ξ_ω, se_ξ_z = cp.se_ξ_z,
          SSR = Optim.minimum(fit), N = length(w))
jldsave(joinpath(@__DIR__, "AggSupply.jld2"); p_star)
println()
println("Saved p_star = ", p_star, " to AggSupply.jld2")

# %% 9. Attach Z, L, λ (already computed in the diagnostics step above, at ρ̂)
# directly onto the estimation panel itself and export under its own name --
# a single self-contained ACS 2001-2024 panel (wages, supply, and task
# aggregates together), not routed through the legacy
# StateAnalysisPreTfp.dta/MakeStateAnalysis.do chain. Not yet wired into
# MakeIRF.do or any other downstream step.
out = DataFrame(statefip = StateAnalysis.statefip, year = StateAnalysis.year,
                 Supply_Domestic = LD, Supply_Foreign = LF,
                 Wage_Domestic = wD, Wage_Foreign = wF,
                 Z = Z, L = L, lambda = λ)
CSV.write(joinpath(data, "StateAggSupplyAcs.csv"), out)
println("Wrote ", nrow(out), " rows to ", joinpath(data, "StateAggSupplyAcs.csv"))
