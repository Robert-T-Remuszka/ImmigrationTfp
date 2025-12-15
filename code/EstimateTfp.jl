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
NS, NT  = length(unique(StateAnalysis[:, :statefip])), length(unique(StateAnalysis[:,:year]))
x0 = vcat(p0.ρ, p0.θ, p0.ζᶠ, Δ0, p0.αᶠ, p0.αᵈ, p0.μ);
#==============================================================================================
VISUALIZATION
1. Are there any restrictions we should place on s̄ for good behavior of the optimization?
==============================================================================================#
#=
How does the cutoff change with the ζ's and the α's? Let's plot some surfaces.
=#
ζ_range = range(0., 20., 100);
Δ_range = range(0., 10., 100);
sbar_vals = [mean(T(AuxParametersConstructor(ζᶠ = ζᶠ, Δ = Δ))[1]) for ζᶠ in ζ_range, Δ in Δ_range];
surface(ζ_range, Δ_range, sbar_vals, xlabel = L"\zeta^F", ylabel = L"\Delta", camera = (35, 25))

#=
How about with the α's?
=#
α_range = range(1e-1, 10., 100);
sbar_vals = [mean(T(AuxParametersConstructor(αᶠ = αᶠ, αᵈ = αᵈ))[1]) for αᶠ in α_range, αᵈ in α_range];
surface(α_range, α_range, sbar_vals)

#=
What does the objective function look like?
=#
ζ_range = range(0., 1., 30);
Δ_range = range(0., 5., 30);
MSE_vals = [SSE(vcat(x0[1:2], [ζᶠ, Δ], x0[5:end])) / (NS * NT) for ζᶠ in ζ_range, Δ in Δ_range];
surface(ζ_range, Δ_range, MSE_vals, xlabel = L"\zeta^F", ylabel = L"\Delta", camera = (30, 30))

α_range = range(1e-1, 3., 50);
MSE_vals = [SSE(vcat(x0[1:4], [αᶠ, αᵈ], x0[end])) / (NS * NT) for αᶠ in α_range, αᵈ in α_range];
surface(α_range, α_range, MSE_vals, label = "MSE", grid = false)

#==============================================================================================
ESTIMATION
==============================================================================================#
res = EstimateProduction(x0);
x_star, MSE_star = res[1], res[2];

# Pack the solution and calculate TFP
p = AuxParametersConstructor(; ρ = x_star[1], θ = x_star[2], ζᶠ = x_star[3], Δ = x_star[4], αᶠ = x_star[5], αᵈ = x_star[6], μ = x_star[end])
StateAnalysis[:, :Z] = Residual(p; df = StateAnalysis)[2]
StateAnalysis[:, :L] = Residual(p; df = StateAnalysis)[3]

# Save
CSV.write(joinpath(data, "StateTfpAndTaskAgg.csv"), StateAnalysis[:,[:statefip, :year, :Z, :L]])