# %% Setup
using DataFrames, StatFiles, TidierData, Statistics, Plots, ForwardDiff, LaTeXStrings, Optim, JLD2, CSV

include("Globals.jl")
include("ProdFunc.jl")

# Load and clean data
StateAnalysis = @chain DataFrame(load(joinpath(data, "StateAnalysisPreTfp.dta"))) begin
    @mutate(
        CapStock = Float64.(CapStock),
        Supply_Foreign = Float64.(Supply_Foreign),
        Supply_Domestic = Float64.(Supply_Domestic),
        GDP = Float64.(GDP),
        Wage_Domestic = Float64.(Wage_Domestic),
        Wage_Foreign = Float64.(Wage_Foreign)
    )
end

# %% Extract arrays for estimation
Y  = StateAnalysis.GDP;
K  = StateAnalysis.CapStock;
wD = StateAnalysis.Wage_Domestic;
wF = StateAnalysis.Wage_Foreign;
LF = StateAnalysis.Supply_Foreign;
LD = StateAnalysis.Supply_Domestic;

w     = wD ./ wF;
max_w = maximum(w);

years   = sort(unique(StateAnalysis.year));
T_years = length(years);
time_id = Int.(indexin(StateAnalysis.year, years));

println("N obs = ", length(Y), "  N years = ", T_years)
println("relative wage w: min=", minimum(w), " max=", maximum(w), " (max_w = ", max_w, ")")

# %% Specification: intercept ιₜ, θ, μ, γ, ρ. (μ,ρ) are bounded as discussed in the paper.
lb = vcat([-Inf, 0.05], [1e-4 / max_w], [0.1, 0.01], fill(-Inf, T_years - 1))
ub = vcat([ Inf, 0.95], [0.999 / max_w], [10.0, 0.49], fill(Inf, T_years - 1))
x0 = vcat([0.0, 0.5], [0.5 / max_w], [1.0, 0.25], zeros(T_years - 1))   # strictly interior to the box; Fminbox requires this
println("x0 length = $(length(x0))")

est = EstimateProdFunc(x0, lb, ub, time_id, T_years, Y, K, wD, wF, LF, LD; show_trace = true)

println()
println("Converged: ", Optim.converged(est.result))
println("Iterations: ", Optim.iterations(est.result))
println("RSS at optimum: ", Optim.minimum(est.result))
println("Fitted: ι=$(est.ι) θ=$(est.θ) γ=$(est.γ) ρ=$(est.ρ)")
println("γ bounds (0.1,10): distance from lower = ", est.γ - 0.1, ", from upper = ", 10.0 - est.γ)
println("ρ bounds (0.01,0.49): distance from lower = ", est.ρ - 0.01, ", from upper = ", 0.49 - est.ρ)
println("μ·max_w: ", est.μ * max_w)
println("ιₜ range: ", extrema(est.ιₜ))

# %% Persist the estimate for downstream use (e.g. Solve_Baseline.jl)
p_star = (ρ = est.ρ, θ = est.θ, γ = est.γ, μ = est.μ, ι = est.ι, ιₜ = est.ιₜ, years = years)
jldsave(joinpath(@__DIR__, "ProductionFunction.jld2"); p_star)
println("Saved p_star = ", p_star, " to ProductionFunction.jld2")

# %% Z variation check
println()
println("=== Z variation ===")
let
    Z_check = getproperty.(TaskAggregates_μ.(est.ρ, est.γ, est.μ, w), :Z)
    println("Z range: ", extrema(Z_check))
    println("Z coefficient of variation: ", std(Z_check) / mean(Z_check))
end

# %% Fitted vs actual, and residual diagnostics
let
    ta = TaskAggregates_μ.(est.ρ, est.γ, est.μ, w)
    Z  = getproperty.(ta, :Z)
    L  = LaborAggregate.(getproperty.(ta, :λ), est.ρ, LF, LD)
    ŷ  = est.ι .+ est.ιₜ[time_id] .+ est.θ .* log.(K) .+ (1 - est.θ) .* (log.(Z) .+ log.(L))

    global fitted_plot = plot(log.(Y), ŷ, seriestype = :scatter, markersize = 2,
        markerstrokewidth = 0, xlabel = "ln Y (actual)", ylabel = "ln Y (fitted)",
        legend = false, title = "Fitted vs Actual")
    plot!(fitted_plot, log.(Y), log.(Y), linewidth = 1, linestyle = :dash, color = :red)
end
fitted_plot

# %% Write Z, L, lambda to StateTfpAndTaskAgg.csv (Handed to MakeStateAnalysis.do -> MakeIRF.do)
let
    ta = TaskAggregates_μ.(est.ρ, est.γ, est.μ, w)
    Z  = getproperty.(ta, :Z)
    λ  = getproperty.(ta, :λ)
    L  = LaborAggregate.(λ, est.ρ, LF, LD)

    out = DataFrame(statefip = StateAnalysis.statefip, year = StateAnalysis.year,
                     Z = Z, L = L, lambda = λ)
    CSV.write(joinpath(data, "StateTfpAndTaskAgg.csv"), out)
    println("Wrote ", nrow(out), " rows to ", joinpath(data, "StateTfpAndTaskAgg.csv"))
end