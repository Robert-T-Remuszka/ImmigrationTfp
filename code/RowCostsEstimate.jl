# %% Setup
using DataFrames, StatFiles, JLD2

include("Globals.jl")
include("RowCostsEstimate_Functions.jl")

# %% Inputs: total 2015 world population by nativity and the five calibration
# targets (two stock shares, three flow rates). total_Lᵈ, total_Lᶠ are
# normalized (total_Lᵈ = 1, total_Lᶠ = the domestic/foreign ratio) rather
# than passed in as raw headcounts -- the steady-state distribution only
# depends on that ratio, not the absolute scale, so this is a free choice,
# and picking a numerically tame one keeps the solver's arithmetic
# well-scaled instead of spanning ~10 orders of magnitude.
tp = load_total_population()
total_Lᵈ = 1.0
total_Lᶠ = tp.total_Lᶠ / tp.total_Lᵈ
(; target_D_abroad, target_D_return, target_F_us, target_F_immigration, target_F_application) = row_cost_targets()

println("total_Lᵈ = ", total_Lᵈ, "  total_Lᶠ = ", total_Lᶠ)
println("target_D_abroad = ", target_D_abroad, "  target_D_return = ", target_D_return)
println("target_F_us = ", target_F_us, "  target_F_immigration = ", target_F_immigration)
println("target_F_application = ", target_F_application)

# %% Admission probability ξ. The realized immigration rate is ξ times the
# application share, so matching both targets pins ξ at their ratio whatever
# the other four parameters are; the remaining four targets then pin
# (fᵈ_ROW, wᵈ_ROW, fᶠ_ROW, wᶠ_ROW) at that ξ.
ξ_gallup = target_F_immigration / target_F_application
println("ξ = ", ξ_gallup)

# %% Calibrate along a path of ξ from 1 (no quota) down to ξ_gallup, each
# solve started from the previous one. The start at ξ = 1 is the home-bias
# calibration's effective Rest of World wages, Wⁿ_ROW^Bⁿ_ROW.
(; Wᵈ_ROW, Wᶠ_ROW) = load_row_income()
rc = load_row_costs()
u = log.([rc.fᵈ_ROW, Wᵈ_ROW^rc.Bᵈ_ROW, rc.fᶠ_ROW, Wᶠ_ROW^rc.Bᶠ_ROW])
W0 = load_wages()
W0ᵈ, W0ᶠ = W0.Wᵈ, W0.Wᶠ

fixed = load_fixed_parameters()
L = fixed.L
ξ_path = [1.0, 0.3, 0.1, 0.03, ξ_gallup]
results = Dict{Float64, NamedTuple}()

for ξ in ξ_path
    println("\n=== Calibrating at ξ = ", ξ, " ===")
    t0 = time()
    cal = calibrate_row(; ξ, total_Lᵈ, total_Lᶠ, target_D_abroad, target_D_return,
                        target_F_us, target_F_immigration, u0 = u, W0ᵈ, W0ᶠ)
    println("done in ", round(time() - t0, digits = 1), "s: ", cal)

    # Verify the fit directly at the calibrated point
    p  = Parameters(fixed, cal.fᵈ_ROW, cal.fᵈ_ROW, cal.fᶠ_ROW, cal.fᶠ_ROW)
    ss = solve_steady_state(p; Wᵈ_ROW = cal.wᵈ_ROW, Wᶠ_ROW = cal.wᶠ_ROW, total_Lᵈ, total_Lᶠ, W0ᵈ, W0ᶠ, ξ)
    m  = row_calibration_moments(ss, p, ξ, total_Lᵈ)
    mean_Wᶠ_us = sum(ss.Lᶠ[1:L] .* ss.Wᶠ[1:L]) / sum(ss.Lᶠ[1:L])
    mean_Vᶠ_us = sum(ss.Lᶠ[1:L] .* ss.Vᶠ[1:L]) / sum(ss.Lᶠ[1:L])
    return_rate_F = sum(ss.Lᶠ[1:L] .* ss.Πᶠ[1:L, p.N]) / sum(ss.Lᶠ[1:L])

    println("--- Fit at ξ = ", ξ, " ---")
    println("share_D_abroad = ", m.share_D_abroad, "  (target ", target_D_abroad, ")")
    println("return_rate_D  = ", m.return_rate_D, "  (target ", target_D_return, ")")
    println("share_F_us     = ", m.share_F_us, "  (target ", target_F_us, ")")
    println("immigration_rate_F = ", m.immigration_rate_F, "  (target ", target_F_immigration, ")")
    println("application_share_F = ", m.application_share_F, "  (target ", target_F_application, ")")
    println("foreign-born return rate = ", return_rate_F)
    println("wᶠ_ROW = ", cal.wᶠ_ROW, "  mean US foreign wage = ", mean_Wᶠ_us,
            "  Vᶠ_ROW - mean Vᶠ_US = ", ss.Vᶠ[p.N] - mean_Vᶠ_us)
    flush(stdout)

    results[ξ] = (; cal..., ξ, target_D_abroad, target_D_return, target_F_us, target_F_immigration,
                    target_F_application,
                    model_D_abroad = m.share_D_abroad, model_D_return = m.return_rate_D,
                    model_F_us = m.share_F_us, model_F_immigration = m.immigration_rate_F,
                    model_F_application = m.application_share_F, return_rate_F, mean_Wᶠ_us,
                    Vᶠ_gap = ss.Vᶠ[p.N] - mean_Vᶠ_us)

    global u = log.([cal.fᵈ_ROW, cal.wᵈ_ROW, cal.fᶠ_ROW, cal.wᶠ_ROW])
    global W0ᵈ, W0ᶠ = ss.Wᵈ[1:L], ss.Wᶠ[1:L]
end

# %% Save the no-quota and Gallup-ξ calibrations, plus the whole ξ path,
# with model-implied moments next to their targets so downstream
# table-building doesn't need to re-solve the model.
jldsave(joinpath(@__DIR__, "RowCostsQuota.jld2");
    baseline = merge(results[ξ_gallup], (; gallup_us_desired)),
    no_quota = results[1.0],
    path     = [results[ξ] for ξ in ξ_path])
println("\nWrote ", joinpath(@__DIR__, "RowCostsQuota.jld2"))
