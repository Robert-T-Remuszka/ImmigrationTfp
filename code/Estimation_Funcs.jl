# estimation of (νᵈ, νᶠ) directly from ACS migration-flow
# data. EstimateScaleBetaAcs.do runs it and saves NuBetaEstimatesAcs.dta. This file's job is packaging
# that Stata/Julia estimation output for Solve_Baseline_Functions.jl's Parameters() to consume
# (load_cp_estimate, load_theta_delta, load_aggsupply_estimate, load_scale_estimate below) --
# it does not re-run any estimation itself.
#
# MakeIRF.do estimates Indirect inference targets and saves the Iv1990 baseline (not the
# Iv1990_LOO robustness check) to IRFEstimates.dta; load_irf_estimates() below
# packages that the same way load_scale_estimate() packages the ν's.

using DataFrames, StatFiles, JLD2

"""
Load (μ_z, ξ_ω, ξ_z) from EstimateCp.do's saved CpEstimates.dta -- the
individual-level choice-probability probit (eq. 4.1). See
notes/ProdFunc Reparamertize v3.md, "Empirical Strategy: Choice Probabilities".
"""
function load_cp_estimate()

    d = only(eachrow(DataFrame(load(joinpath(data, "CpEstimates.dta")))))

    return (μ_z = d.mu_z_hat, ξ_ω = d.xi_w_hat, ξ_z = d.xi_z_hat,
            se_μ_z = d.mu_z_se, se_ξ_ω = d.xi_w_se, se_ξ_z = d.xi_z_se, N = d.N)

end

"""
Load the state-level capital share θ_l and depreciation rate δ_l from
CalibrateTheta.do's saved ThetaDelta.dta. θ_l is backed out of the capital
FOC (3.9), not estimated -- see that file's header for the exact construction
(1994-2021 time average of El-Shagi-Yamarik's state-year deprate).
"""
load_theta_delta() = DataFrame(load(joinpath(data, "ThetaDelta.dta")))

"""
Load the aggregate-supply parameter estimate (ρ, plus the (μ_z,ξ_ω,ξ_z) it was
fit taking as given) from AggSupply_Estimate.jl's saved AggSupply.jld2.
"""
load_aggsupply_estimate() = load(joinpath(@__DIR__, "AggSupply.jld2"), "p_star")

"""
Load (νᵈ, νᶠ) from EstimateScaleBetaAcs.do's saved NuBetaEstimatesAcs.dta.
Uses the raw-wage OLS row for each nativity: that file's Durbin-Wu-Hausman
test fails to reject OLS==IV on either equation, and OLS is far more precise
for no detectable loss of consistency, so OLS is the "directly inferred"
estimate reported in the paper (see EstimateScaleBetaAcs.do's header for
the full diagnosis of why 2SLS is unusable at this identification strategy's
power).
"""
function load_scale_estimate()

    df = DataFrame(load(joinpath(data, "NuBetaEstimatesAcs.dta")))
    d  = only(eachrow(df[df.nativity .== "D", :]))
    f  = only(eachrow(df[df.nativity .== "F", :]))

    return (νᵈ = d.nu_ols, νᶠ = f.nu_ols)

end

"""
Load the section-4 LPIV impulse responses (Iv1990 baseline only) from
MakeIRF.do's saved IRFEstimates.dta. Returns a Dict keyed by outcome name
("Z", "L", "Wage_Domestic", "Wage_Foreign"), each holding the horizon-ordered
vectors (h, β, se, F) needed to target these moments in indirect inference.
"""
function load_irf_estimates()

    df = DataFrame(load(joinpath(data, "IRFEstimates.dta")))
    sort!(df, [:outcome, :h])

    return Dict(
        first(g.outcome) => (h = g.h, β = g.beta, se = g.se, F = g.Fstat)
        for g in groupby(df, :outcome)
    )

end
