using DataFrames, StatFiles, LinearAlgebra, Plots, LaTeXStrings, CSV, TidierData, Statistics, Optim, ForwardDiff
include("Globals.jl");
include("Functions.jl");

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

p0 = AuxParametersConstructor();
Δ0 = 5.;
δ0 = 1.;
T = length(unique(StateAnalysis[:, :year]));
N = length(unique(StateAnalysis[:, :statefip]));
x0 = vcat(p0.ρ, p0.θ, p0.ζᶠ, Δ0, δ0, p0.ξ, p0.χ, p0.Inter);
#==============================================================================================
VISUALIZATION
1. Are there any restrictions we should place on s̄ for good behavior of the optimization?
==============================================================================================#
#=
How does s̄ change with the ζ's and the α's? Let's plot some surfaces.
=#
ζ_range = range(0., 20., 100);
Δ_range = range(0., 10., 100);
sbar_vals = [mean(s̄(AuxParametersConstructor(ζᶠ = ζᶠ, Δ = Δ))[1]) for ζᶠ in ζ_range, Δ in Δ_range];
surface(ζ_range, Δ_range, sbar_vals, xlabel = L"\zeta^F", ylabel = L"\Delta", camera = (35, 25))

δ_range = range(-2., 10., 50); # REMARK: Below 2 and s̄ > 1 which leads to numerical instability so we will avoid that.
sbar_vals = [mean(s̄(AuxParametersConstructor(δ = δ))[1]) for δ in δ_range];
scatter(δ_range, sbar_vals, xlabel = L"\delta", grid = false, label = L"\bar s")

#=
What does the objective function looks like?
=#
ζ_range = range(0., 5., 30);
Δ_range = range(0., 5., 30);
MSE_vals = [SSE(vcat(x0[1:2], [ζᶠ, Δ], x0[5:end])) / (N * T) for ζᶠ in ζ_range, Δ in Δ_range];
surface(ζ_range, Δ_range, MSE_vals, xlabel = L"\zeta^F", ylabel = L"\Delta", camera = (70, 30))

δ_range = range(-2., 10., 50);
MSE_vals = [SSE(vcat(x0[1:4], [δ], x0[6:end])) / (N * T) for δ in δ_range];
scatter(δ_range, MSE_vals, xlabel = L"\delta", label = "MSE", grid = false)

#==============================================================================================
ESTIMATION
==============================================================================================#
res = EstimateProduction(x0);
x_star, MSE_star = res[1], res[2];

# Pack the solution and calculate TFP
p = AuxParametersConstructor(;ρ = x_star[1], θ = x_star[2], ζᶠ = x_star[3], Δ = x_star[4], δ = x_star[5], ξ = x_star[6: 4 + T], χ = x_star[5 + T : 3 + T + N], Inter = x_star[end])
StateAnalysis[:, :Z] = Residual(p; df = StateAnalysis)[2]
StateAnalysis[:, :L] = Residual(p; df = StateAnalysis)[3]

# Save
CSV.write(joinpath(data, "StateTfpAndTaskAgg.csv"), StateAnalysis[:,[:statefip, :year, :Z, :L]])