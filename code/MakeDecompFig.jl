# Illustrative theory figure for Section 2.2 (imperfect specialization):
# producers draw idiosyncratic comparative advantage, so instead of a sharp
# task boundary, the foreign-born task share λ(w) responds smoothly to the
# relative wage w := W^D/W^F. AggSupply_Functions.jl's TaskAggregates_LN
# implements the paper's closed form for this (Appendix B.2) directly, so
# this script calls it rather than re-deriving anything:
#   TaskAggregates_LN(ρ, μ_z, ξ_ω, ξ_z, w) -> (; Z, λ, oneMinusλ)
# where Z is aggregate task productivity and λ the foreign-born task share.
# `oneMinusλ` is returned separately from λ because naively computing 1-λ
# loses precision once λ saturates near 0 or 1 -- see that function's own
# docstring.
#
# This script illustrates three facts about that object:
#   (1) The market-clearing relative wage w is increasing in relative foreign
#       labor supply s := L^F/L^D -- equivalently, the foreign/domestic wage
#       ratio W^F/W^D = 1/w is *decreasing* in s: more relative foreign labor
#       supply pushes foreign-born wages down relative to domestic. Market
#       clearing requires
#         w^{1/(1-ρ)} = s · (1-λ(w)) / λ(w),
#       whose two sides are monotonic in w (LHS increasing, RHS decreasing
#       since λ(w) is increasing), so it has a unique solution, found here by
#       bisection.
#   (2) Z(w) is hump-shaped in the relative wage w, with a unique interior
#       maximizer (found here numerically -- it need not be exactly w=1 in
#       general, only when ξ_ω, ξ_z, μ_z happen to balance that way). The
#       calibrated steady state's state-level wage ratios (SteadyState.jld2)
#       cluster almost exactly there.
#   (3) The foreign-born task share λ(w) itself is smoothly increasing in w --
#       the mechanism underlying both (1) and (2): a higher relative domestic
#       wage shifts tasks toward foreign-born labor.
#
# Parameters are the actual calibrated values -- ρ from AggSupply_Estimate.jl's
# AggSupply.jld2, (μ_z, ξ_ω, ξ_z) from EstimateCp.do's CpEstimates.dta --
# loaded via the same Estimation_Funcs.jl helpers MakeTables.jl uses, not
# retyped from the paper's printed table. Axes are left numerically unlabeled
# throughout: this is a conceptual illustration of *shape*, not a
# calibrated-magnitude plot (the one exception is the shaded band on the Z(w)
# panel, which marks the actual calibrated steady-state wage-ratio range from
# SteadyState.jld2).

using CairoMakie, StatsFuns, JLD2, Optim

include("Globals.jl")
include("Estimation_Funcs.jl")
include("AggSupply_Functions.jl")

# %% Load the calibrated parameters
aggsupply = load_aggsupply_estimate()
cp        = load_cp_estimate()

ρ   = Float64(aggsupply.ρ)
μ_z = Float64(cp.μ_z)
ξ_ω = Float64(cp.ξ_ω)
ξ_z = Float64(cp.ξ_z)

println("Calibrated parameters used: ρ = $ρ, μ_z = $μ_z, ξ_ω = $ξ_ω, ξ_z = $ξ_z")

# %% Closed-form objects (paper's Appendix B.2 identity, via AggSupply_Functions.jl)
Z_of(w)     = TaskAggregates_LN(ρ, μ_z, ξ_ω, ξ_z, w).Z
λ_of(w)     = TaskAggregates_LN(ρ, μ_z, ξ_ω, ξ_z, w).λ
oneMλ_of(w) = TaskAggregates_LN(ρ, μ_z, ξ_ω, ξ_z, w).oneMinusλ

# Productivity-maximizing relative wage, found numerically (bounded Brent
# search -- Z(w) is unimodal, so this is well-posed).
opt   = optimize(w -> -Z_of(w), 0.1, 5.0)
wstar = Optim.minimizer(opt)
Zstar = -Optim.minimum(opt)

# Market-clearing relative wage w(s), solved by bisection on (1e-4, 20): both
# sides of w^{1/(1-ρ)} = s·(1-λ(w))/λ(w) are monotonic in w, so this is
# well-posed without a dedicated root-finding package.
function solve_w(s; lo = 1e-4, hi = 20.0, tol = 1e-12, maxit = 200)
    f(w) = w^(1 / (1 - ρ)) - s * oneMλ_of(w) / λ_of(w)
    a, c = lo, hi
    fa = f(a)
    for _ in 1:maxit
        (c - a) < tol && break
        m = (a + c) / 2
        fm = f(m)
        if sign(fm) == sign(fa)
            a, fa = m, fm
        else
            c = m
        end
    end
    return (a + c) / 2
end

# Illustrative relative-supply grid, s = L^F/L^D ≈ 0.05 to 0.6. The current US
# calibration (foreign-born ≈18.7% of the labor force) is s ≈ 0.187/0.813.
s_current     = 0.187 / (1 - 0.187)
s_low, s_high = 0.05, 0.6
s_grid        = range(s_low, s_high, length = 200)
w_grid        = solve_w.(s_grid)

s_markers = (s_low, s_current, s_high)
w_markers = solve_w.(s_markers)

# Calibrated steady-state wage-ratio range (SteadyState.jld2's state-level
# W^D/W^F), for the shaded band on the Z(w) panel. One state (DC) is an
# extreme outlier (w≈0.44) that is not representative of the bulk
# distribution, so the band reported is the range of the other 51 states.
ss_wage_ratios = let d = JLD2.load(joinpath(@__DIR__, "SteadyState.jld2"))
    sort(d["ss"].Wᵈ ./ d["ss"].Wᶠ)
end
ss_band = (ss_wage_ratios[2], ss_wage_ratios[end])

println("w* (Z-maximizing relative wage) = $wstar, Z(w*) = $Zstar")
for (s, w) in zip(s_markers, w_markers)
    println("  s = $(round(s, digits = 3)) -> w = $(round(w, digits = 4)), 1/w = $(round(1/w, digits = 4)), Z(w) = $(round(Z_of(w), digits = 4))")
end
println("Calibrated steady-state wage-ratio range (excl. 1 outlier state): ", ss_band)

# %% Figure
try
    set_theme!(fontsize = 16,
        fonts = (; regular = "Nimbus Roman", bold = "Nimbus Roman", italic = "Nimbus Roman"))
catch
    set_theme!(fontsize = 16)
end

navy = colorant"#1d3557"
gray = colorant"#8d99ae"
band = colorant"#a8dadc"
marker_colors = cgrad(:YlOrRd, 4, categorical = true)[2:4]   # low -> mid -> high, sequential warm palette
marker_labels = ["Low Lᶠ/Lᴰ", "Current Lᶠ/Lᴰ", "High Lᶠ/Lᴰ"]

wlo, whi = 0.3, 3.0
wgrid    = exp.(range(log(wlo), log(whi), length = 400))

fig = Figure(size = (900, 700))

# --- Row 1, left: W^F/W^D = 1/w(s), decreasing in s -------------------------
ax1 = Axis(fig[1, 1], xlabel = "Lᶠ/Lᴰ", ylabel = "Wᶠ/Wᴰ = 1/w")
lines!(ax1, s_grid, 1 ./ w_grid, color = navy, linewidth = 2.5)
for (s, w, c) in zip(s_markers, w_markers, marker_colors)
    scatter!(ax1, [s], [1 / w], color = c, markersize = 12, strokecolor = :black, strokewidth = 0.5)
end

# --- Row 1, right: Z(w) hump, peak marked, calibrated band -----------------
ax2 = Axis(fig[1, 2], xlabel = "w = Wᴰ/Wᶠ", ylabel = "Z")
vspan!(ax2, ss_band[1], ss_band[2], color = (band, 0.4))
lines!(ax2, wgrid, Z_of.(wgrid), color = navy, linewidth = 2.5)
vlines!(ax2, [wstar], color = gray, linestyle = :dash, linewidth = 1.5)
scatter!(ax2, [wstar], [Zstar], color = gray, markersize = 12, strokecolor = :black, strokewidth = 0.5)
for (w, c) in zip(w_markers, marker_colors)
    scatter!(ax2, [w], [Z_of(w)], color = c, markersize = 10, strokecolor = :black, strokewidth = 0.5)
end

# --- Row 2: λ(w), the foreign-born task share -------------------------------
ax3 = Axis(fig[2, 1:2], xlabel = "w = Wᴰ/Wᶠ", ylabel = "λ")
lines!(ax3, wgrid, λ_of.(wgrid), color = navy, linewidth = 2.5)
for (w, c) in zip(w_markers, marker_colors)
    scatter!(ax3, [w], [λ_of(w)], color = c, markersize = 12, strokecolor = :black, strokewidth = 0.5)
end

for ax in (ax1, ax2, ax3)
    hidedecorations!(ax, label = false, ticklabels = true, ticks = true, minorticks = true)
    hidespines!(ax, :t, :r)
end

# Shared legend: the shaded band and the w=1 reference line (both from panel
# 2) alongside the three s-scenario markers used across all three panels.
legend_elements = [
    PolyElement(color = (band, 0.4)),
    LineElement(color = gray, linestyle = :dash, linewidth = 1.5),
    [MarkerElement(color = c, marker = :circle, markersize = 12, strokecolor = :black, strokewidth = 0.5) for c in marker_colors]...,
]
legend_labels = ["Observed real-wage variation in sample", "w = 1", marker_labels...]

Legend(fig[3, 1:2], legend_elements, legend_labels;
    orientation = :horizontal, nbanks = 1, tellwidth = false, tellheight = true, framevisible = false)

outfile = joinpath(graphs, "TaskSpecializationIllustration.pdf")
save(outfile, fig)
println("Saved figure to $outfile ($(filesize(outfile)) bytes)")
