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

# %% Calibration steps. The first has no quota (ι = Inf), so the application
# share is untargeted. The quota steps then lower the application-share
# target to the Gallup one through the shares at which ξ = 0.3, 0.1 and 0.03
# would match the immigration-rate target, each solve started from the
# previous one. The first start is the home-bias calibration's effective
# Rest of World wages, Wⁿ_ROW^Bⁿ_ROW.
steps = vcat([("no quota", nothing)],
             [("ξ ≈ $ξ", target_F_immigration / ξ) for ξ in (0.3, 0.1, 0.03)],
             [("Gallup", target_F_application)])

(; Wᵈ_ROW, Wᶠ_ROW) = load_row_income()
rc = load_row_costs()
u = log.([rc.fᵈ_ROW, Wᵈ_ROW^rc.Bᵈ_ROW, rc.fᶠ_ROW, Wᶠ_ROW^rc.Bᶠ_ROW])
W0 = load_wages()
W0ᵈ, W0ᶠ = W0.Wᵈ, W0.Wᶠ
ξ0 = 1.0

fixed = load_fixed_parameters()
L = fixed.L
results = NamedTuple[]

for (label, target_app) in steps
    println("\n=== Calibrating: ", label, ", application-share target = ", target_app, " ===")

    # The first quota step starts ι 1% below the no-quota inflow, so the quota
    # binds from the start (at the no-quota inflow ξ sits on the kink of
    # min{1, ι/applicants}, where raising ι has no effect), and ξ at the
    # ratio of the two flow targets.
    if target_app !== nothing && length(u) == 4
        global u  = vcat(u, log(0.99 * results[end].ι_inflow))
        global ξ0 = target_F_immigration / target_app
    end

    t0 = time()
    cal = calibrate_row(; total_Lᵈ, total_Lᶠ, target_D_abroad, target_D_return, target_F_us,
                        target_F_immigration, target_F_application = target_app, u0 = u, W0ᵈ, W0ᶠ, ξ0)
    println("done in ", round(time() - t0, digits = 1), "s: ", cal)

    # Verify the fit directly at the calibrated point
    p  = Parameters(fixed, cal.fᵈ_ROW, cal.fᵈ_ROW, cal.fᶠ_ROW, cal.fᶠ_ROW)
    ss = solve_steady_state(p; Wᵈ_ROW = cal.wᵈ_ROW, Wᶠ_ROW = cal.wᶠ_ROW, total_Lᵈ, total_Lᶠ, W0ᵈ, W0ᶠ,
                            ι = cal.ι, ξ0)
    m  = row_calibration_moments(ss, p, total_Lᵈ)
    ι_inflow = ss.ξ * applicants(ss, p)
    mean_Wᶠ_us = sum(ss.Lᶠ[1:L] .* ss.Wᶠ[1:L]) / sum(ss.Lᶠ[1:L])
    mean_Vᶠ_us = sum(ss.Lᶠ[1:L] .* ss.Vᶠ[1:L]) / sum(ss.Lᶠ[1:L])
    return_rate_F = sum(ss.Lᶠ[1:L] .* ss.Πᶠ[1:L, p.N]) / sum(ss.Lᶠ[1:L])

    println("--- Fit: ", label, " ---")
    println("share_D_abroad = ", m.share_D_abroad, "  (target ", target_D_abroad, ")")
    println("return_rate_D  = ", m.return_rate_D, "  (target ", target_D_return, ")")
    println("share_F_us     = ", m.share_F_us, "  (target ", target_F_us, ")")
    println("immigration_rate_F = ", m.immigration_rate_F, "  (target ", target_F_immigration, ")")
    println("application_share_F = ", m.application_share_F, "  (target ", target_app, ")")
    println("ι = ", cal.ι, "  realized inflow = ", ι_inflow, "  (", ι_inflow * tp.total_Lᵈ, " people)  ξ = ", ss.ξ)
    println("foreign-born return rate = ", return_rate_F)
    println("wᶠ_ROW = ", cal.wᶠ_ROW, "  mean US foreign wage = ", mean_Wᶠ_us,
            "  Vᶠ_ROW - mean Vᶠ_US = ", ss.Vᶠ[p.N] - mean_Vᶠ_us)
    flush(stdout)

    push!(results, (; label, cal..., ξ = ss.ξ, ι_inflow, ι_inflow_people = ι_inflow * tp.total_Lᵈ,
                    target_D_abroad, target_D_return, target_F_us, target_F_immigration,
                    target_F_application = target_app,
                    model_D_abroad = m.share_D_abroad, model_D_return = m.return_rate_D,
                    model_F_us = m.share_F_us, model_F_immigration = m.immigration_rate_F,
                    model_F_application = m.application_share_F, return_rate_F, mean_Wᶠ_us,
                    Vᶠ_gap = ss.Vᶠ[p.N] - mean_Vᶠ_us))

    global u  = target_app === nothing ? log.([cal.fᵈ_ROW, cal.wᵈ_ROW, cal.fᶠ_ROW, cal.wᶠ_ROW]) :
                                         log.([cal.fᵈ_ROW, cal.wᵈ_ROW, cal.fᶠ_ROW, cal.wᶠ_ROW, cal.ι])
    global W0ᵈ, W0ᶠ = ss.Wᵈ[1:L], ss.Wᶠ[1:L]
    global ξ0 = ss.ξ
end

# %% Save the no-quota and Gallup calibrations, plus every step, with
# model-implied moments next to their targets so downstream table-building
# doesn't need to re-solve the model. ι is in the model's population units
# (total_Lᵈ = 1); ι_inflow_people converts the realized inflow to people.
jldsave(joinpath(@__DIR__, "RowCostsQuota.jld2");
    baseline = merge(results[end], (; gallup_us_desired)),
    no_quota = results[1],
    path     = results)
println("\nWrote ", joinpath(@__DIR__, "RowCostsQuota.jld2"))
