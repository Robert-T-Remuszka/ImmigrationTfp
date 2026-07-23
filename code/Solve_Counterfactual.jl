# %% Setup
using JLD2, StatFiles, DataFrames, NonlinearSolve, LinearAlgebra, ForwardDiff, Plots

include("Globals.jl");
include("ProdFunc.jl");
include("Solve_Baseline_Functions.jl");
include("Solve_Counterfactual_Functions.jl");

Init_Data = load_init_data();
p = Parameters(; Init_Data);

# %% Solve the baseline first — every counterfactual is defined relative to it
Baseline = SolveBaseline(p);

# %% Solve the counterfactual: households learn at t=1 of a one-standard-deviation MIT
# shock to the US-specific mobility cost (σ=1 normalization), decaying per the AR(1) law
# of motion thereafter
M̂ = [mit_shock(σ = 1.0, ψ = p.ψ)(t) for t in 1:(Baseline.T - 1)];
CF = SolveCounterfactual(M̂; Baseline, p);

# %% Plot the model-implied impulse response: employment-weighted domestic/foreign wage
# and TFP (Z) in the counterfactual relative to the baseline, for US locations
us = 1:(p.N - 1)
horizon = 2:min(42, Baseline.T)   # period 1 is the pre-shock anchor, identical to Baseline; the shock first hits period 2
t = 0:(length(horizon) - 1)       # relabel period 2 as t=0, since that's when the shock first hits

wᵈ_hat = vec(sum(CF.Wᵈ[us, :] .* Baseline.Lᵈ[us, :], dims = 1) ./ sum(Baseline.Lᵈ[us, :], dims = 1)) ./
         vec(sum(Baseline.Wᵈ[us, :] .* Baseline.Lᵈ[us, :], dims = 1) ./ sum(Baseline.Lᵈ[us, :], dims = 1))
wᶠ_hat = vec(sum(CF.Wᶠ[us, :] .* Baseline.Lᶠ[us, :], dims = 1) ./ sum(Baseline.Lᶠ[us, :], dims = 1)) ./
         vec(sum(Baseline.Wᶠ[us, :] .* Baseline.Lᶠ[us, :], dims = 1) ./ sum(Baseline.Lᶠ[us, :], dims = 1))

Z_Baseline = ComputeZ(Baseline; p)
Z_CF       = ComputeZ(CF; p)
Z_hat      = vec(sum(Z_CF[us, :] .* Baseline.Lᵈ[us, :], dims = 1) ./ sum(Baseline.Lᵈ[us, :], dims = 1)) ./
             vec(sum(Z_Baseline[us, :] .* Baseline.Lᵈ[us, :], dims = 1) ./ sum(Baseline.Lᵈ[us, :], dims = 1))

share_hat = vec(sum(CF.Lᶠ[us, :], dims = 1) ./ sum(CF.Lᵈ[us, :] .+ CF.Lᶠ[us, :], dims = 1)) ./
            vec(sum(Baseline.Lᶠ[us, :], dims = 1) ./ sum(Baseline.Lᵈ[us, :] .+ Baseline.Lᶠ[us, :], dims = 1))

p1 = plot(t, wᵈ_hat[horizon], xlabel = "Years since shock", ylabel = "ŵᴰ", title = "Domestic Wage Response", legend = false,
linewidth = 2.5, grid = false);
p2 = plot(t, wᶠ_hat[horizon], xlabel = "Years since shock", ylabel = "ŵᶠ", title = "Foreign Wage Response", legend = false,
linewidth = 2.5, grid = false);
p3 = plot(t, Z_hat[horizon], xlabel = "Years since shock", ylabel = "Ẑ", title = "Productivity Response", legend = false,
linewidth = 2.5, grid = false);
p4 = plot(t, share_hat[horizon], xlabel = "Years since shock", ylabel = "L̂ᶠ/(Lᴰ+Lᶠ)", title = "Foreign-Born Share Response", legend = false,
linewidth = 2.5, grid = false);

irf_plot = plot(p1, p2, p3, p4, layout = (2, 2), size = (1000, 700))
savefig(irf_plot, joinpath(graphs, "CounterfactualIRF.pdf"))
irf_plot