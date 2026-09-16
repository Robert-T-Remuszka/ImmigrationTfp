# %% Setup
using JLD2, StatFiles, DataFrames, NonlinearSolve, LinearAlgebra, ForwardDiff, Plots

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

panels = map((:Z, :L, :Wᵈ, :Wᶠ)) do outcome
    pl = plot(title = String(outcome), xlabel = "Horizon k", ylabel = "Kernel(k)", legend = :outertopright, grid = false)
    for σ_p in sort(collect(keys(kernels)))
        plot!(pl, 0:size(kernels[σ_p][outcome], 2) - 1, aggregate(kernels[σ_p][outcome]), label = "σ_p=$σ_p", linewidth = 1.5)
    end
    pl
end
scalability_plot = plot(panels..., layout = (2, 2), size = (1200, 800))
savefig(scalability_plot, joinpath(graphs, "KernelScalability.pdf"))
scalability_plot
