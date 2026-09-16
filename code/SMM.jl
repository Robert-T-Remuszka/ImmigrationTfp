# %% Setup
using JLD2, StatFiles, DataFrames, NonlinearSolve, LinearAlgebra, ForwardDiff, Plots, Random, Statistics

include("Globals.jl")
include("Estimation_Funcs.jl")
include("AggSupply_Functions.jl")
include("Solve_Baseline_Functions.jl")
include("Solve_Counterfactual_Functions.jl")
include("SMM_Functions.jl")

Init_Data = load_init_data()
p = Parameters(; Init_Data)

# %% Solve for the model's own steady state, re-anchored away from the 1996 transition --
# see SolveSteadyState's docstring for why shocking the raw 1996 path directly would fail
# BKM's time-invariance requirement
ss = SolveSteadyState(p)
println("Steady state reached at t★ = $(ss.t★) (transition T = $(ss.Transition.T)); ",
        "re-anchored Baseline_ss flatness = $(ss.flatness) after $(ss.iters) refinement(s)")
jldsave(joinpath(@__DIR__, "SMM_SteadyState.jld2");
        p_ss = ss.p_ss, Baseline_ss = ss.Baseline_ss, t★ = ss.t★, flatness = ss.flatness)

# %% BKM scalability/sign-antisymmetry check: does the normalized kernel 𝒦(σ_p) = IRF(σ_p)/σ_p
# stay put as the probe size σ_p grows, and is it antisymmetric in sign? This is BKM's own
# accuracy check (their §5.2, Figs. 5-6), and it is what licenses treating any single probe's
# IRF as "the" linear kernel reused for every future shock draw.
kernels, report = kernel_scalability(ss.Baseline_ss; p = ss.p_ss, K = 60)
show(report, allrows = true)

# %% Plot: normalized kernel overlay per outcome, all probes -- the direct analogue of BKM's
# Fig. 5/6. Aggregated (employment-weighted, at the steady-state level) across US locations
# for legibility.
us = 1:ss.p_ss.N - 1
weights = ss.Baseline_ss.Lᵈ[us, 1] .+ ss.Baseline_ss.Lᶠ[us, 1]
aggregate(κ_x) = vec(sum(κ_x .* weights, dims = 1) ./ sum(weights))

panels = map((:Lᵈ, :Lᶠ, :Wᵈ, :Wᶠ, :Z, :L)) do outcome
    pl = plot(title = String(outcome), xlabel = "Horizon k", ylabel = "Kernel(k)", legend = :outertopright, grid = false)
    for σ_p in sort(collect(keys(kernels)))
        plot!(pl, 0:size(kernels[σ_p][outcome], 2) - 1, aggregate(kernels[σ_p][outcome]), label = "σ_p=$σ_p", linewidth = 1.5)
    end
    pl
end
scalability_plot = plot(panels..., layout = (3, 2), size = (1200, 1200))
savefig(scalability_plot, joinpath(graphs, "KernelScalability.pdf"))
scalability_plot

# %% Build the σ_p=0.01 kernel (the chosen probe, per the scalability check above) and smoke-test
# the BKM convolution: draw a shared national innovation sequence once, then simulate a panel at
# a placeholder σ (the actual grid search over (σ,ψ) is the next step, not this one).
κ, _ = bkm_kernel(ss.Baseline_ss; p = ss.p_ss, σ_p = 0.01, K = 60)

Tsim = 30
e    = draw_innovations(Tsim; seed = 1)
σ_placeholder = 0.05
panel = simulate_panel(κ, e; σ = σ_placeholder)

println("Simulated panel shapes: ", NamedTuple(k => size(v) for (k, v) in pairs(panel)))
println("L: std across (state,time) = ", round(std(panel.L), digits = 4),
        "  range = ", round.(extrema(panel.L), digits = 4))

# %% Plot: simulated L path for a handful of states, to eyeball that the convolution produces
# sensible-looking dynamics (persistent, decaying responses to the shared shock, not white noise).
sim_plot = plot(title = "Simulated log L deviation (σ=$σ_placeholder)", xlabel = "Simulated period",
                ylabel = "log L̂", legend = :outertopright, grid = false)
for l in 1:5:size(panel.L, 1)
    plot!(sim_plot, 1:Tsim, panel.L[l, :], label = "state $l", linewidth = 1.5)
end
sim_plot
