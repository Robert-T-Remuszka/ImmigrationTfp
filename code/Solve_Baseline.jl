# %% Setup
using JLD2, StatFiles, DataFrames, NonlinearSolve, LinearAlgebra, ForwardDiff, Plots

include("Globals.jl")
include("ProdFunc.jl")
include("Solve_Baseline_Functions.jl")

Init_Data = load_init_data()
p = Parameters(; Init_Data)

# %% Solve the baseline transition to steady state, holding the US-specific mobility cost
# mₜ constant at its long-run mean forever (no_shock — see Solve_Baseline_Functions.jl)
Baseline = SolveBaseline(p)

# %% Plot employment-weighted transition paths from the 1996 initial data
t = 0:(Baseline.T - 1)                 # years since 1996 (t = 0 is the 1996 GDP-anchored period)
us = 1:p.N - 1                     # the 50 states + DC, excluding the exogenous "Rest of World" row

wᵈ_agg = vec(sum(Baseline.Wᵈ[us, :] .* Baseline.Lᵈ[us, :], dims = 1) ./ sum(Baseline.Lᵈ[us, :], dims = 1))
wᶠ_agg = vec(sum(Baseline.Wᶠ[us, :] .* Baseline.Lᶠ[us, :], dims = 1) ./ sum(Baseline.Lᶠ[us, :], dims = 1))

share_Dᵈ_us = vec(sum(Baseline.Lᵈ[us, :], dims = 1) ./ sum(Baseline.Lᵈ, dims = 1)) .* 100   # % of domestic-born (US citizens) residing in the US
share_Lᶠ_us = vec(sum(Baseline.Lᶠ[us, :], dims = 1) ./ sum(Baseline.Lᶠ, dims = 1)) .* 100   # % of foreign-born residing in the US (immigrant stock share)

p1 = plot(t, wᵈ_agg, xlabel = "Years since 1996", ylabel = "\$ per worker", title = "Employment-Weighted Domestic Wage", legend = false,
linewidth = 1.5, grid = false);
p2 = plot(t, wᶠ_agg, xlabel = "Years since 1996", ylabel = "\$ per worker", title = "Employment-Weighted Foreign Wage", legend = false,
linewidth = 1.5, grid = false);
p3 = plot(t, share_Dᵈ_us, xlabel = "Years since 1996", ylabel = "%", title = "Domestic-Born Population Share in US", legend = false,
linewidth = 1.5, grid = false);
p4 = plot(t, share_Lᶠ_us, xlabel = "Years since 1996", ylabel = "%", title = "Foreign-Born Population Share in US", legend = false,
linewidth = 1.5, grid = false);

transition_plot = plot(p1, p2, p3, p4, layout = (2, 2), size = (1000, 700), linewidth = 2.5)
savefig(transition_plot, joinpath(graphs, "BaselineTransition.pdf"))
transition_plot