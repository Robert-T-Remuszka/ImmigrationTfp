# %% Setup
using DataFrames, StatFiles, JLD2

include("Globals.jl")
include("RowCostsEstimate_Functions.jl")

# %% Inputs: real *measured* ROW wages (before any home-bias adjustment),
# total 2015 world population by nativity, and the four calibration targets
# (two stock shares, two flow rates). total_Lᵈ, total_Lᶠ are normalized
# (total_Lᵈ = 1, total_Lᶠ = the domestic/foreign ratio) rather than passed in
# as raw headcounts -- the steady-state distribution only depends on that
# ratio, not the absolute scale, so this is a free choice, and picking a
# numerically tame one keeps the solver's arithmetic well-scaled instead of
# spanning ~10 orders of magnitude.
(; Wᵈ_ROW, Wᶠ_ROW) = load_row_income()
Wᵈ_ROW_meas, Wᶠ_ROW_meas = Wᵈ_ROW, Wᶠ_ROW
tp = load_total_population()
total_Lᵈ = 1.0
total_Lᶠ = tp.total_Lᶠ / tp.total_Lᵈ
(; target_D_abroad, target_D_return, target_F_us, target_F_immigration) = row_cost_targets()

println("Wᵈ_ROW_meas = ", Wᵈ_ROW_meas, "  Wᶠ_ROW_meas = ", Wᶠ_ROW_meas)
println("total_Lᵈ = ", total_Lᵈ, "  total_Lᶠ = ", total_Lᶠ)
println("target_D_abroad = ", target_D_abroad, "  target_D_return = ", target_D_return)
println("target_F_us = ", target_F_us, "  target_F_immigration = ", target_F_immigration)

# %% Calibrate the two symmetric ROW costs and two home-bias amenities
# jointly. Seeded near a region a manual diagnostic search already confirmed
# is smooth and well-behaved (see project notes) -- no grid search needed
# here, unlike the abandoned asymmetric-cost/fixed-wage approach, which
# never had a region like this to seed from in the first place.
println("Calibrating (fᵈ_ROW, Bᵈ_ROW, fᶠ_ROW, Bᶠ_ROW)...")
t0 = time()
(; fᵈ_ROW, Bᵈ_ROW, fᶠ_ROW, Bᶠ_ROW) = calibrate_home_bias(;
    Wᵈ_ROW_meas, Wᶠ_ROW_meas, total_Lᵈ, total_Lᶠ,
    target_D_abroad, target_D_return, target_F_us, target_F_immigration)
println("Calibration done in ", round(time() - t0, digits = 1), "s")

println("fᵈ_ROW = ", fᵈ_ROW, "  Bᵈ_ROW = ", Bᵈ_ROW)
println("fᶠ_ROW = ", fᶠ_ROW, "  Bᶠ_ROW = ", Bᶠ_ROW)

# %% Verify the fit directly at the calibrated point
fixed = load_fixed_parameters()
p  = Parameters(fixed, fᵈ_ROW, fᵈ_ROW, fᶠ_ROW, fᶠ_ROW)
Wᵈ_ROW_eff = Wᵈ_ROW_meas^Bᵈ_ROW
Wᶠ_ROW_eff = Wᶠ_ROW_meas^Bᶠ_ROW
ss = solve_steady_state(p; Wᵈ_ROW = Wᵈ_ROW_eff, Wᶠ_ROW = Wᶠ_ROW_eff, total_Lᵈ, total_Lᶠ,
    W0ᵈ = load_wages().Wᵈ, W0ᶠ = load_wages().Wᶠ)

share_D_abroad = ss.Lᵈ[p.N] / total_Lᵈ
return_rate_D  = leave_probability(ss.Vᵈ, p.fᵈ, p.β, p.νᵈ, p.N)
Lᵈ_us_ss = sum(ss.Lᵈ[1:p.N - 1])
Lᶠ_us_ss = sum(ss.Lᶠ[1:p.N - 1])
share_F_us         = Lᶠ_us_ss / (Lᵈ_us_ss + Lᶠ_us_ss)
immigration_rate_F = leave_probability(ss.Vᶠ, p.fᶠ, p.β, p.νᶠ, p.N)

println("--- Fit at calibrated point ---")
println("share_D_abroad = ", share_D_abroad, "  (target ", target_D_abroad, ")")
println("return_rate_D  = ", return_rate_D, "  (target ", target_D_return, ")")
println("share_F_us     = ", share_F_us, "  (target ", target_F_us, ")")
println("immigration_rate_F = ", immigration_rate_F, "  (target ", target_F_immigration, ")")

# %% Save alongside the other estimated-parameter jld2 files (task shares,
# scale parameters, ρ) -- model-implied moments are saved alongside their
# targets so downstream table-building doesn't need to re-solve the model.
jldsave(joinpath(@__DIR__, "RowCosts.jld2");
    p_star = (; fᵈ_ROW, Bᵈ_ROW, fᶠ_ROW, Bᶠ_ROW,
                target_D_abroad, target_D_return, target_F_us, target_F_immigration,
                model_D_abroad = share_D_abroad, model_D_return = return_rate_D,
                model_F_us = share_F_us, model_F_immigration = immigration_rate_F))
