using JLD2, DataFrames, TidierData, Plots, LaTeXStrings
include("ProdFunc.jl");
include("PF_Transition.jl");
@load "ProductionFunction.jld2" p_star;

T = 20;
Params1 = Parameters(p_star; T = T);
m  = [-Params1.ψ^t for t in 0:T-2];        # Impulse response of AR(1) to one-unit reudcition in US migration costs
μ̇  = exp.(-(1 - Params1.ψ) .* m);          # Proportional changes in exp(m_t)
Soln = SolveTransition(μ̇; p = Params1);

# Check if u̇ goes to one
plot(1:T, Soln.U̇ᵈ[1,:], grid = false, linewidth = 3., 
label = L"$\dot u^d_{t+1}$", dpi = 300, xticks = 1:T, title = "Lifetime Utilility Changes")
plot!(1:T, Soln.U̇ᶠ[1,:], grid = false, linewidth = 3.,label = L"$\dot u^f_{t+1}$")