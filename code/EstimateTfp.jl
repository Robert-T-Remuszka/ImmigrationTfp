using DataFrames, StatFiles, TidierData, LinearAlgebra, Plots, LaTeXStrings, CSV
using Optimization, Zygote
include("Globals.jl");
include("Functions.jl");

# Read in the state level data
StateAnalysis = DataFrame(load(joinpath(data, "StateAnalysisPreTfp.dta")))

p0 = AuxParameters();
T = length(unique(StateAnalysis[:, :year]));
N = length(unique(StateAnalysis[:,:statefip]));
x0 = vcat(p0.ρ, p0.θ, p0.αᶠ, p0.αᵈ, p0.Inter, p0.δ[2:end], 2. * p0.ξ[2:end], 1. * p0.ζᶠ, 6. * p0.ζᵈ);

#==============================================================================================
VISUALIZATION
Vary each parameter one by one from the starting point to see how the residuals change.
==============================================================================================#

#=
For each of the scalar parameters, loop through a grid and plot
=#
plots = [];
labels = [L"\alpha^F", L"\alpha^D", L"\zeta^F", L"\zeta^D"];
params = [3, 4, length(x0) - 1, length(x0)];
for (param, label) in zip(params, labels)
    
    upper = x0[param] + 5.;
    lower = x0[param] - .5;
    
    vals = range(lower, upper, length = 10);
    yvals = []
    for v in vals
        x = copy(x0)
        x[param] = v
        push!(yvals, SSE(x; df = StateAnalysis) / (N * T))
    end

    push!(plots, scatter(vals, yvals, xlabel = label, grid = false, linewidth = 2., legend = false, ylabel = "MSE"))

end

plot(plots...)

#==============================================================================================
ESTIMATION
==============================================================================================#
g(x, p) = SSE(x)
f = OptimizationFunction(g, Optimization.AutoForwardDiff())

# Set bounds
lb = zeros(5 + N - 1 + T - 1 + 2);               # Initialize all bounds first
lb[5] = -Inf;                                    # Inter
lb[6: 5 + N - 1] .= -Inf;                        # SFEs
lb[5 + N : 5 + N - 1 + T - 1] .= -Inf;           # TFEs

ub = ones(5 + N - 1 + T - 1 + 2);                 # Initialize all bounds first
ub[3] = Inf;                                      # αᶠ
ub[4] = Inf;                                      # αᵈ
ub[5] = Inf;                                      # Inter
ub[6: 5 + N - 1] .= Inf;                          # SFEs
ub[5 + N : 5 + N - 1 + T - 1] .= Inf;             # TFEs
ub[end-1:end] .= Inf;                             # Comp adv params

prob = OptimizationProblem(f, x0, lb = lb, ub = ub);
sol = solve(prob, Optimization.LBFGS())

# Pack the solution and calculate TFP
p = AuxParameters(sol[1], sol[2], sol[3], sol[4], sol[5], vcat(0., sol[6: 5 + N - 1]), vcat(0., sol[5 + N : 5 + N - 1 + T - 1]), sol[end - 1], sol[end]);
StateAnalysis[:, :Z] = Residual(p; df = StateAnalysis)[2]
StateAnalysis[:, :L] = Residual(p; df = StateAnalysis)[3]

# Save
CSV.write(joinpath(data, "StateTfpAndTaskAgg.csv"), StateAnalysis[:,[:statefip, :year, :Z, :L]])