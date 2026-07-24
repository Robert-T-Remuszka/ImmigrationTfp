# %% Setup
using JLD2, StatFiles, DataFrames, CSV, NonlinearSolve, LinearAlgebra, ForwardDiff, Plots, Random, Optim

include("Globals.jl")
include("ProdFunc.jl")
include("Solve_Baseline_Functions.jl")
include("Solve_Counterfactual_Functions.jl")
include("Simulation_Functions.jl")
include("IndInf_Functions.jl")

# %% Load the state panel underlying the paper's LPIV (built by Functions.do's
# PreRegProcessing), forcing state/year to Int for clean joins in add_shift.
RegData = DataFrame(load(joinpath(data, "StateAnalysisRegReady.dta")))
RegData.state = Int.(RegData.state)
RegData.year  = Int.(RegData.year)

# %% Validate: run Julia LPIV/2SLS against the actual data and compare to the
# Stata (ivreg2) coefficients exported by MakeIRF.do, for all four target variables.
depvars = [:Z, :Wage_Domestic, :Wage_Foreign, :L]

for v in depvars

    β_julia = lpiv(RegData, v)
    Stata   = CSV.read(joinpath(data, "BaselineIRF_$(v).csv"), DataFrame)
    β_stata = Stata.Beta_Iv1990

    maxdiff = maximum(abs.(β_julia .- β_stata))
    println("$v: max |Julia - Stata| = $maxdiff")

end
# %% Load the empirical target (coefficients + SEs, for inverse-variance weighting) that
# indirect inference will match
Init_Data = load_init_data()

data_β  = (; (v => CSV.read(joinpath(data, "BaselineIRF_$(v).csv"), DataFrame).Beta_Iv1990 for v in depvars)...)
data_se = (; (v => CSV.read(joinpath(data, "BaselineIRF_$(v).csv"), DataFrame).Se_Iv1990   for v in depvars)...)

# %% Validate the full simulated-moments pipeline at a single (default) parameter vector
# before attempting any optimization — one Baseline+counterfactual solve, S small draws
seeds = collect(1:10)
sim_β_test = simulated_moments(4.5, 4.5, 0.50, 1.0; Init_Data, seeds, verbose = false)

for v in depvars
    println("$v: data β[1:3] = ", round.(data_β[v][1:3], digits = 4),
            "   sim β[1:3] = ", round.(sim_β_test[v][1:3], digits = 4))
end

Q_test = smd_objective(4.5, 4.5, 0.50, 1.0; Init_Data, data_β, data_se, seeds)
println("SMD objective at default θ: ", Q_test)

# %% Full indirect-inference estimation: search over (νᴰ, νᶠ, ψ, σ) to minimize the
# inverse-variance-weighted SMD objective against the real Figure-1 LPIV coefficients.
# Uses S common-random-number draws (fixed across all θ) so the objective is smooth/
# deterministic in θ, per Boppart, Krusell, and Mitman (2018). Each evaluation solves the
# baseline + one counterfactual from scratch (~2 minutes), so this is bounded by a
# wall-clock time limit rather than an iteration count, returning the best θ found so far.
S_estimation = 20
seeds_estimation = collect(1:S_estimation)

θ0    = [4.5, 4.5, 0.50, 1.0]     # νᴰ, νᶠ, ψ, σ starting guess
lower = [0.5, 0.5, -0.95, 0.01]
upper = [20.0, 20.0, 0.95, 5.0]

function obj(θ)
    νᴰ, νᶠ, ψ, σ = θ
    Q = smd_objective(νᴰ, νᶠ, ψ, σ; Init_Data, data_β, data_se, seeds = seeds_estimation)
    println("θ = $θ   Q = $Q")
    flush(stdout)
    return Q
end

result = Optim.optimize(obj, lower, upper, θ0, Fminbox(NelderMead()),
    Optim.Options(time_limit = 3 * 3600.0, show_trace = false))

θ̂ = Optim.minimizer(result)
νᴰ̂, νᶠ̂, ψ̂, σ̂ = θ̂
println("Estimated θ̂: νᴰ=$νᴰ̂, νᶠ=$νᶠ̂, ψ=$ψ̂, σ=$σ̂   (Q=$(Optim.minimum(result)))")

jldsave(joinpath(data, "IndInfEstimate.jld2"); θ̂, νᴰ̂, νᶠ̂, ψ̂, σ̂, Q = Optim.minimum(result))

# %% Plot the estimated model's IRF against the real Figure-1 LPIV coefficients (with the
# empirical ±1.645·SE band) — the key output of this estimation.
sim_β_hat = simulated_moments(νᴰ̂, νᶠ̂, ψ̂, σ̂; Init_Data, seeds = seeds_estimation)

ylabs = Dict(:Z => "Z", :Wage_Domestic => "wᴰ", :Wage_Foreign => "wᶠ", :L => "L")
plots = map(depvars) do v
    h = 0:9
    lo = data_β[v] .- 1.645 .* data_se[v]
    hi = data_β[v] .+ 1.645 .* data_se[v]
    plt = plot(h, data_β[v], ribbon = (data_β[v] .- lo, hi .- data_β[v]),
        fillalpha = 0.2, label = "Data", linewidth = 2, color = :steelblue,
        xlabel = "h", ylabel = "Δʰln($(ylabs[v]))", title = string(v), grid = false)
    plot!(plt, h, sim_β_hat[v], label = "Model", linewidth = 2, color = :firebrick, linestyle = :dash)
    plt
end

irf_comparison = plot(plots..., layout = (2, 2), size = (1000, 700))
savefig(irf_comparison, joinpath(graphs, "ModelVsDataIRF.pdf"))
irf_comparison
