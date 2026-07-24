# %% Setup
using JLD2, StatFiles, DataFrames, NonlinearSolve, LinearAlgebra, ForwardDiff, Plots, Random

include("Globals.jl");
include("ProdFunc.jl");
include("Solve_Baseline_Functions.jl");
include("Solve_Counterfactual_Functions.jl");
include("Simulation_Functions.jl");

Init_Data = load_init_data();
p = Parameters(; Init_Data);

# %% Solve the baseline, then compute each US location's own unit-shock IRF.
Baseline = SolveBaseline(p);
IRF = compute_state_irf(Baseline; p, K = 60);

# %% Simulate a panel: draw one shared national innovation sequence and convolve each
# state's own IRF with it (BKM).
σ = 1.0
Tsim = 30
sim, e = simulate_panel(IRF, σ, Tsim; seed = 1234)

# %% Showcase: plot a few example states' simulated domestic-wage paths against the shared
# innovation sequence driving them, to see the convolution mechanism at work
us = 1:(p.N - 1)
K = size(IRF.wᴰ, 2)
states_to_plot = [1, 5, 20]

p1 = plot(1:Tsim, sim.wᴰ[states_to_plot, :]', xlabel = "Simulated period", ylabel = "log ŵᴰ",
    title = "Simulated Domestic Wage Response (example states)",
    label = permutedims(string.("State ", states_to_plot)), linewidth = 1.5, grid = false)
p2 = plot(1:Tsim, e[K:(Tsim + K - 1)], xlabel = "Simulated period", ylabel = "eₜ",
    title = "Driving Innovation Sequence (national)", legend = false, linewidth = 1, grid = false)

showcase_plot = plot(p1, p2, layout = (2, 1), size = (900, 700))
savefig(showcase_plot, joinpath(graphs, "SimulationShowcase.pdf"))
showcase_plot