# CDP (2019) eq.-(27)-style estimation of (νᵈ, νᶠ) directly from ACS migration-flow
# data. Nu_Derivation.md derives the estimating equation; EstimateScaleBetaAcs.do
# runs it (pair-FE IV, wage ratio instrumented with its own lag and a lag-2
# flow-ratio term) and saves NuBetaEstimatesAcs.dta. This file's job is packaging
# that Stata output for Julia consumption, mirroring how Solve_Baseline_Functions.jl's
# load_prodfunc_estimate() packages ProdFunc_Estimate.jl's output -- it does not
# re-run any estimation itself.
#
# ψ (the AR(1) persistence of the US-specific mobility cost mₜ) and σ (its
# innovation SD) are not estimated this way: they're identified only through
# their effect on the dynamics of migration flows, wages and productivity, so
# they require indirect inference against the paper's section-4 LPIV impulse
# responses. That Stata estimation doesn't exist yet -- load_irf_estimates()
# below is a placeholder for whatever eventually produces it.

using DataFrames, StatFiles

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
Placeholder for the future handoff of the section-4 LPIV impulse responses
(Z, L, wᴰ, wᶠ to a migration shock) into Julia, for indirect-inference
estimation of (ψ, σ). Not yet built -- no Stata script produces these IRFs
for this paper yet (distinct from the older MakeIRF.do, which estimates a
different paper's Bartik/TFP design).
"""
function load_irf_estimates()

    error("load_irf_estimates: not yet implemented -- the section 4 LPIV Stata script doesn't exist yet")

end
