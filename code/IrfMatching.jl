# %% Setup
using DataFrames, StatFiles, JLD2, Plots

include("Globals.jl")
include("IrfMatching_Functions.jl")

default(dpi = 300)

# %% Build the calibrated Parameters (symmetric ROW costs + home-bias
# amenities from RowCostsEstimate.jl) and solve the pre-shock steady state.
fixed = load_fixed_parameters()
(; fᵈ_ROW, Bᵈ_ROW, fᶠ_ROW, Bᶠ_ROW) = load_row_costs()
p = Parameters(fixed, fᵈ_ROW, fᵈ_ROW, fᶠ_ROW, fᶠ_ROW)

(; Wᵈ_ROW, Wᶠ_ROW) = load_row_income()
Wᵈ_ROW_eff = Wᵈ_ROW^Bᵈ_ROW
Wᶠ_ROW_eff = Wᶠ_ROW^Bᶠ_ROW

tp = load_total_population()
total_Lᵈ = 1.0
total_Lᶠ = tp.total_Lᶠ / tp.total_Lᵈ

ss = solve_steady_state(p; Wᵈ_ROW = Wᵈ_ROW_eff, Wᶠ_ROW = Wᶠ_ROW_eff, total_Lᵈ, total_Lᶠ,
    W0ᵈ = load_wages().Wᵈ, W0ᶠ = load_wages().Wᶠ)

jldsave(joinpath(@__DIR__, "SteadyState.jld2"); ss, p, total_Lᵈ, total_Lᶠ)
println("Steady state solved and saved to SteadyState.jld2")

# %% Solve the shooting problem for a placeholder (ψ, σ): a one-standard-
# deviation *reduction* in the foreign-born ROW→US migration cost (σ<0),
# decaying at rate ψ. These are NOT yet the calibrated values -- ψ, σ are
# what this whole exercise will estimate; the point here is only to confirm
# the shooting algorithm itself works before calibrating anything. Sized as
# 10% of the calibrated fᶠ_ROW, purely so the shock is a plausible order of
# magnitude for a first test.
ψ0 = 0.5
σ0 = -0.1 * fᶠ_ROW
T0 = 300

println("Solving shooting problem: ψ=$ψ0, σ=$σ0, T=$T0...")
t0 = time()
path = solve_transition(p, ss; ψ = ψ0, σ = σ0, T = T0)
println("Done in ", round(time() - t0, digits = 1), "s")

# %% Sanity checks: did the path actually return to the steady state by T?
Lᵈ_gap = maximum(abs.(path.Lᵈ[:, end] .- ss.Lᵈ) ./ ss.Lᵈ)
Lᶠ_gap = maximum(abs.(path.Lᶠ[:, end] .- ss.Lᶠ) ./ ss.Lᶠ)
println("Max relative gap to steady state at t=T: Lᵈ=", Lᵈ_gap, "  Lᶠ=", Lᶠ_gap)
if Lᵈ_gap > 1e-4 || Lᶠ_gap > 1e-4
    println("WARNING: gap is large -- T may not be long enough, consider increasing T0")
end

# %% Model-implied IRFs for the same four outcomes MakeIRF.do's empirical
# LPIV targets (Z, L, Wage_Foreign, Wage_Domestic): log-deviation of the
# shooting path's national aggregate from the steady state at each horizon --
# this is the object to compare against the empirical β_h (long-differenced
# specification), which Piger & Stockwell (2025) show estimates the same
# level-relative-to-no-shock-counterfactual object, just more efficiently in
# small samples than a levels specification would.
irf = national_irfs(path, ss, p)
h   = 0:T0

irf_plot = plot(
    plot(h, irf.Z,  lw = 2, legend = false, xlabel = "h", ylabel = "Δ log Z",  title = "Z"),
    plot(h, irf.L,  lw = 2, legend = false, xlabel = "h", ylabel = "Δ log L",  title = "L"),
    plot(h, irf.WF, lw = 2, legend = false, xlabel = "h", ylabel = "Δ log wᶠ", title = "Wage Foreign"),
    plot(h, irf.WD, lw = 2, legend = false, xlabel = "h", ylabel = "Δ log wᴰ", title = "Wage Domestic"),
    layout = (2, 2), size = (1000, 700),
    plot_title = "Model-Implied IRFs: One-SD Reduction in Foreign ROW→US Migration Cost"
)
savefig(irf_plot, joinpath(graphs, "IrfMatching_ModelIRFs.pdf"))
irf_plot
