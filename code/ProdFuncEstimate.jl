using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using DataFrames, StatFiles, LinearAlgebra, LaTeXStrings, CSV, TidierData, Statistics, Optim, ForwardDiff
using Plots, JLD2

include("Globals.jl");
include("ProdFunc.jl");

# Read in the state level data
StateAnalysis = @chain DataFrame(load(joinpath(data, "StateAnalysisPreTfp.dta"))) begin
    @mutate(
        CapStock = Float64.(CapStock),
        Supply_Foreign = Float64.(Supply_Foreign), 
        Supply_Domestic = Float64.(Supply_Domestic),
        GDP = Float64.(GDP),
        Wage_Domestic = Float64.(Wage_Domestic),
        Wage_Foreign = Float64.(Wage_Foreign)
    )
end;

p0 = AuxParameters();
N, T = length(unique(StateAnalysis[:, :statefip])), length(unique(StateAnalysis[:,:year]));
x0 = vcat(p0.ρ, p0.θ, p0.γᶠ, p0.Δ, p0.Γ, p0.αᶠ, p0.αᵈ, p0.ψ, p0.ι, p0.ιₛ[2 : end], p0.ιₜ[2: end]); # Fixed effects have leading zeros

#=======
VISUALIZATION
=======#

#= RESIDUAL SUM OF SQUARES 

Now look at how the residual sum of squares changes with the parameters of the 
production function.
=#

# ι
ι_range  = range(-10., 10., 100);
plot(ι_range, ι -> RSS(vcat(x0[1:8], ι, x0[10:end])), grid = false, linewidth = 2., xlabel = L"\iota", legend = false)

# (γᶠ, Δ)
γᶠ_range = range(1e-3, 2., 20);
Δ_range  = range(1e-1, 2., 20);
heatmap(γᶠ_range, Δ_range, (γᶠ, Δ) -> RSS(vcat(x0[1:2], γᶠ, Δ, x0[5:end])), grid = false, linewidth = 2., xlabel = L"\Delta", ylabel = L"\gamma^F", c=:viridis)

# Γ
Γ_range = range(-5, 5, 100);
plot(Γ_range, Γ -> RSS(vcat(x0[1:4], Γ, x0[6:end])), grid = false, linewidth = 2., xlabel = L"\Gamma", legend = false)

# ρ
ρ_range = range(0.1, 0.99, 100);
plot(ρ_range, ρ -> RSS(vcat(ρ, x0[2:end])), grid = false, linewidth = 2., xlabel = L"\rho", legend = false)

# θ
θ_range = range(0., 1., 100);
plot(θ_range, θ -> RSS(vcat(x0[1], θ, x0[3:end])), grid = false, linewidth = 2., xlabel = L"\theta", legend = false)

# (αᶠ, αᵈ)
α_range = range(1e-4, 10., 20);
heatmap(α_range, α_range, (αᶠ, αᵈ) -> RSS(vcat(x0[1:5], αᶠ, αᵈ, x0[8:end])), xlabel = L"\alpha^F", ylabel = L"\alpha^D")

#=======
ESTIMATION
=======#
x_star, MSE_star = EstimateProdFunc(x0);
p_star  = AuxParameters(ρ = x_star[1], θ = x_star[2], γᶠ = x_star[3], Δ = x_star[4], Γ = x_star[5], αᶠ = x_star[6], αᵈ = x_star[7], 
ψ = x_star[8], ι = x_star[9], ιₛ = vcat(0., x_star[10 : 8 + N]), ιₜ = vcat(0., x_star[9 + N : end]));
@save "ProductionFunction.jld2" p_star;

# Add estimated objects from production function to data
ComputeZ(wᶠ, wᵈ, F, D) = ComputeReduced(wᶠ, wᵈ, F, D; p = p_star).Z;
ComputeL(wᶠ, wᵈ, F, D) = ComputeReduced(wᶠ, wᵈ, F, D; p = p_star).L;
Computeλ(wᶠ, wᵈ, F, D) = ComputeReduced(wᶠ, wᵈ, F, D; p = p_star).λ;
StateAnalysis = @chain StateAnalysis begin
   @mutate(
    Z = ComputeZ(Wage_Foreign, Wage_Domestic, Supply_Foreign, Supply_Domestic),
    L = ComputeL(Wage_Foreign, Wage_Domestic, Supply_Foreign, Supply_Domestic),
    lambda = Computeλ(Wage_Foreign, Wage_Domestic, Supply_Foreign, Supply_Domestic)
   )  
end

# Save the estimates
CSV.write(joinpath(data, "StateTfpAndTaskAgg.csv"), StateAnalysis[:,[:statefip, :year, :Z, :L, :lambda]]);