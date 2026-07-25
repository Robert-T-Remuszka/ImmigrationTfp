using JLD2, StatFiles, DataFrames, CSV, NonlinearSolve, LinearAlgebra, ForwardDiff, Random, Plots

include(joinpath(@__DIR__, "Globals.jl"))
include(joinpath(@__DIR__, "ProdFunc.jl"))
include(joinpath(@__DIR__, "Solve_Baseline_Functions.jl"))
include(joinpath(@__DIR__, "Solve_Counterfactual_Functions.jl"))
include(joinpath(@__DIR__, "Simulation_Functions.jl"))
include(joinpath(@__DIR__, "IndInf_Functions.jl"))

# %% Load the empirical target and the currently-saved parameter estimate (whatever θ̂ is
# on disk right now — a single serial evaluation, no need for the Distributed workers used
# by the outer search in IndInf_Estimate.jl).
Init_Data = load_init_data()
depvars = [:Z, :Wage_Domestic, :Wage_Foreign, :L]

data_β  = (; (v => CSV.read(joinpath(data, "BaselineIRF_$(v).csv"), DataFrame).Beta_Iv1990 for v in depvars)...)
data_se = (; (v => CSV.read(joinpath(data, "BaselineIRF_$(v).csv"), DataFrame).Se_Iv1990   for v in depvars)...)

# Placeholder parameters (not the final estimate) — the best point found so far by the
# parallel global search of the still-running indirect-inference job (3865225), targeting
# only the Z and L IRFs with equal weighting; local refinement (xNES) is still ongoing.
νᴰ̂, νᶠ̂, ψ̂, σ̂ = 1.4724, 1.2582, -0.8704, 1.8112
println("Plotting at global-search-best placeholder θ (indirect inference still pending): νᴰ=$νᴰ̂, νᶠ=$νᶠ̂, ψ=$ψ̂, σ=$σ̂")
flush(stdout)

S_estimation = 20
seeds_estimation = collect(1:S_estimation)

sim_β_hat = simulated_moments(νᴰ̂, νᶠ̂, ψ̂, σ̂; Init_Data, seeds = seeds_estimation)

# %% Plot the model's IRF against the real Figure-1 LPIV coefficients (±1.645·SE band).
# Z, L (the targeted moments) share the top row; the two untargeted wage IRFs share the
# bottom row. Model line is orange throughout — solid where targeted, dashed where not — so
# the distinction reads at a glance without needing separate colors.
plot_order = [:Z, :L, :Wage_Domestic, :Wage_Foreign]
ylabs = Dict(:Z => "Z", :Wage_Domestic => "wᴰ", :Wage_Foreign => "wᶠ", :L => "L")
plots = map(plot_order) do v
    h = 0:9
    lo = data_β[v] .- 1.645 .* data_se[v]
    hi = data_β[v] .+ 1.645 .* data_se[v]
    targeted = v in TARGET_DEPVARS
    plt = plot(h, data_β[v], ribbon = (data_β[v] .- lo, hi .- data_β[v]),
        fillalpha = 0.2, label = "Data", linewidth = 2, color = :steelblue,
        xlabel = "h", ylabel = "Δʰln($(ylabs[v]))", title = string(v), grid = false,
        left_margin = 8Plots.mm)
    plot!(plt, h, sim_β_hat[v], label = targeted ? "Model (targeted)" : "Model (untargeted)",
        linewidth = 2, color = :orange, linestyle = targeted ? :solid : :dash)
    plt
end

irf_comparison = plot(plots..., layout = (2, 2), size = (1000, 700))
savefig(irf_comparison, joinpath(graphs, "ModelVsDataIRF.pdf"))
println("Wrote plot to $(joinpath(graphs, "ModelVsDataIRF.pdf"))")