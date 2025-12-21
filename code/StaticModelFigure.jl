using DataFrames, StatFiles, LaTeXStrings, TidierData, Statistics
using CairoMakie
using Colors: RGB

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

p0 = AuxParameters(ρ = 0.50, θ = 0.50, γᶠ = 1., Δ = 2., Γ = -1., αᶠ = 1., αᵈ = 2., ι = 9.);

#=======
VISUALIZATION
=======#
# Take a look at the comparative advantage schedule - show how productivity effects depend on parameter Γ
z(τ)  = exp(p0.Γ + (p0.Δ) * τ);
z⁰(τ) = exp((p0.Δ) * τ);
τ_range = range(0, 1, 100);

fig = Figure();
ax  = Axis(fig[1,1], xlabel = L"\tau", limits = (0, 1, 0, 3));
hlines!([1.], color = :black, alpha = 0.5);
lines!(ax, τ_range, z.(τ_range),  linewidth = 2., color = RGB(0/255, 147/255, 245/255), label = L"\Gamma < 0");
lines!(ax, τ_range, z⁰.(τ_range), linewidth = 2., color = RGB(247/255, 129/255, 4/255), label = L"\Gamma = 0");
hideydecorations!(ax, ticklabels = false, ticks = false);
hidexdecorations!(ax, label = false, ticks = false, ticklabels = false);
hidespines!(:t, :r);
axislegend(ax, position = :rb);
fig

#= PRODUCTIVITY V. MARGINAL TASK

Now look at productivty. To do this with the functions we will use to estimate the production function,
normalize wᶠ = 1 so that the relative wage is given by wᵈ. Moving wᵈ around will move 𝒯 around and this
is what we want to plot Z against.

=#
mean_rel_wage = mean(StateAnalysis[:, :Wage_Domestic] ./ StateAnalysis[:, :Wage_Foreign]);
wage_range    = range(0., mean_rel_wage + 2., 100);
Z(wᵈ) = ComputeReduced(1., wᵈ, 1., 1.; p = p0).Z;
τ(wᵈ) = ComputeReduced(1., wᵈ, 1., 1.; p = p0).𝒯;

fig = Figure();
ax1 = Axis(fig[1,1], xlabel = L"\mathcal{T}");
ax2 = Axis(fig[1,1], yaxisposition = :right);
Prod = lines!(ax1, τ.(wage_range), Z.(wage_range), linewidth = 2., color = RGB(0/255, 147/255, 245/255));
CA   = lines!(ax2, τ.(wage_range), z.(τ.(wage_range)), linewidth = 2., color = RGB(247/255, 129/255, 4/255));
hlines!([1.], color = :black, alpha = 0.5, linestyle = :dash);
vlines!([-p0.Γ/(p0.Δ)], color = :black, alpha = 0.5, linestyle = :dash);
hideydecorations!(ax1);
hideydecorations!(ax2);
hidexdecorations!(ax1, label = false);
hidexdecorations!(ax2, label = false);
hidespines!(ax1, :t, :r);
hidespines!(ax2, :t, :r);
Legend(fig[1, 1], [Prod, CA], [L"Z(\mathcal{T})", L"z^D(\mathcal{T})/z^F(\mathcal{T})"]; 
tellheight = false, tellwidth = false, halign = :left, valign = :top, margin = (10, 10, 10, 10), framevisible = true)
fig
save("../output/graphs/ProductivityAndMarginalTasks.pdf", fig);
