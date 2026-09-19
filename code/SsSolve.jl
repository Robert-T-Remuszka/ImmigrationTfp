# %% Setup
using DataFrames, StatFiles, JLD2, Plots, Random

include("Globals.jl");
include("SsSolve_Functions.jl");

default(dpi = 300);

# %% fᵈ_ROW, fᶠ_ROW are not yet calibrated -- NaN placeholders until then.
p = Parameters(; fᵈ_ROW = NaN, fᶠ_ROW = NaN);

# %% Partial-equilibrium demo: solve Vᵈ, Vᶠ, Πᵈ, Πᶠ at a placeholder wage
# vector, with fᵈ_ROW, fᶠ_ROW set to 0. Wages aren't solved for yet, so a
# little noise is added around 1 just to see some variation in this exercise
# -- not a real wage vector.
Random.seed!(1);
p_pe = Parameters(; fᵈ_ROW = 0.0, fᶠ_ROW = 0.0);
Wᵈ = exp.(0.3 .* randn(p_pe.N));
Wᶠ = exp.(0.3 .* randn(p_pe.N));
s = Solution(p_pe, Wᵈ, Wᶠ);

# State name labels (abbreviated), read off StateAnalysisPreTfp.dta and sorted
# to match load_bilateral_costs()'s FIPS-code state ordering.
names_df = unique(select(DataFrame(load(joinpath(data, "StateAnalysisPreTfp.dta"))), [:statefip, :state]));
sort!(names_df, :statefip);
labels = [[get(state_abbrevs, n, n) for n in names_df.state]; "ROW"];

# %% Value function by location. y-axis is zoomed to the actual range of V --
# on a from-zero axis the cross-location variation is imperceptible next to
# V's overall level.
value_plot = plot(
    bar(1:p_pe.N, s.Vᵈ, title = "Vᵈ", legend = false, xticks = (1:p_pe.N, labels),
        xrotation = 45, xtickfontsize = 5, grid = false,
        ylims = extrema(s.Vᵈ) .+ (-0.1, 0.1) .* (extrema(s.Vᵈ)[2] - extrema(s.Vᵈ)[1])),
    bar(1:p_pe.N, s.Vᶠ, title = "Vᶠ", legend = false, xticks = (1:p_pe.N, labels),
        xrotation = 45, xtickfontsize = 5, grid = false,
        ylims = extrema(s.Vᶠ) .+ (-0.1, 0.1) .* (extrema(s.Vᶠ)[2] - extrema(s.Vᶠ)[1])),
    layout = (2, 1), size = (1000, 700)
)

# %% Choice probabilities: full Πᵈ, Πᶠ as heatmaps (log scale -- the diagonal
# stay probability is orders of magnitude larger than any single off-diagonal
# entry, so a linear color scale hides the cross-state pattern entirely)
choiceprob_heatmap = plot(
    heatmap(log10.(s.Πᵈ), title = "log₁₀ Πᵈ", xlabel = "Destination", ylabel = "Origin", yflip = true, grid = false),
    heatmap(log10.(s.Πᶠ), title = "log₁₀ Πᶠ", xlabel = "Destination", ylabel = "Origin", yflip = true, grid = false),
    layout = (1, 2), size = (1100, 450)
)

# %% Heatmaps of the inverted bilateral migration costs fᵈ, fᶠ, for the paper.
# Both panels share one color scale so the domestic/foreign cost structures
# are directly comparable.
costs = load_bilateral_costs()
state_labels = labels[1:end - 1]
clims = (0, max(maximum(costs.fᵈ), maximum(costs.fᶠ)))

# Tight axis limits (cells are centered on 1:51) so there's no default padding
# between the axis and the first/last row or column of the heatmap.
tight_lims = (0.5, 51.5)

d_heatmap = heatmap(costs.fᵈ, title = "Domestic Born", xlabel = "Destination", ylabel = "Origin",
    xticks = (1:51, state_labels), yticks = (1:51, state_labels), xrotation = 90,
    xtickfontsize = 7, ytickfontsize = 7, yflip = true, clims = clims,
    xlims = tight_lims, ylims = tight_lims,
    color = :viridis, aspect_ratio = 1, colorbar = false, grid = false,
    left_margin = 12Plots.mm, bottom_margin = 8Plots.mm)
f_heatmap = heatmap(costs.fᶠ, title = "Foreign Born", xlabel = "Destination", ylabel = "Origin",
    xticks = (1:51, state_labels), yticks = (1:51, state_labels), xrotation = 90,
    xtickfontsize = 7, ytickfontsize = 7, yflip = true, clims = clims,
    xlims = tight_lims, ylims = tight_lims,
    color = :viridis, aspect_ratio = 1, colorbar = false, grid = false,
    left_margin = 12Plots.mm, bottom_margin = 8Plots.mm, right_margin = 10Plots.mm)

# A single shared colorbar in its own panel, so it doesn't shrink either
# heatmap's plotting area the way putting colorbar=true on one panel would.
# All-NaN heatmap draws no visible cells, only the colorbar itself.
cbar = heatmap(fill(NaN, 2, 2), clims = clims, color = :viridis, colorbar = true,
    colorbar_title = "Cost", framestyle = :none, legend = false,
    xaxis = false, yaxis = false, xticks = false, yticks = false,
    right_margin = 10Plots.mm)

cost_heatmap = plot(d_heatmap, f_heatmap, cbar,
    layout = @layout([a b c{0.10w}]), size = (2000, 850))
savefig(cost_heatmap, joinpath(graphs, "BilateralCostsHeatmap.pdf"))
cost_heatmap
