# %% Setup
using DataFrames, StatFiles, JLD2

include("Globals.jl")
include("RowCostsEstimate_Functions.jl")

#=================================================================
Autarky counterfactual, steady state to steady state. The baseline is the
calibrated non-stochastic steady state. The counterfactual makes entering the
US infinitely costly for the foreign-born: row N of fᶠ (Rest of World as
origin) is Inf for every US destination and 0 on the diagonal. The
foreign-born already in the US still draw taste shocks and leave, and none
can enter, so the US has no foreign-born labor in the new steady state.
=================================================================#

# %% Baseline: calibrated parameters and steady state
fixed = load_fixed_parameters()
(; fᵈ_ROW, Bᵈ_ROW, fᶠ_ROW, Bᶠ_ROW) = load_row_costs()
p = Parameters(fixed, fᵈ_ROW, fᵈ_ROW, fᶠ_ROW, fᶠ_ROW)

# Home bias enters through the effective Rest of World wage, Wⁿ_ROW^Bⁿ_ROW
(; Wᵈ_ROW, Wᶠ_ROW) = load_row_income()
Wᵈ_ROW_eff = Wᵈ_ROW^Bᵈ_ROW
Wᶠ_ROW_eff = Wᶠ_ROW^Bᶠ_ROW

tp = load_total_population()
total_Lᵈ = 1.0
total_Lᶠ = tp.total_Lᶠ / tp.total_Lᵈ

t0 = time()
ss = solve_steady_state(p; Wᵈ_ROW = Wᵈ_ROW_eff, Wᶠ_ROW = Wᶠ_ROW_eff, total_Lᵈ, total_Lᶠ,
    W0ᵈ = load_wages().Wᵈ, W0ᶠ = load_wages().Wᶠ)
println("Baseline steady state solved in ", round(time() - t0, digits = 1), "s")

L = p.N - 1
println("Baseline foreign-born share of US labor = ", sum(ss.Lᶠ[1:L]) / (sum(ss.Lᶠ[1:L]) + sum(ss.Lᵈ[1:L])))

# %% Autarky: entering the US is infinitely costly for the foreign-born
fᶠ_autarky = copy(p.fᶠ)
fᶠ_autarky[p.N, 1:L] .= Inf
fᶠ_autarky[p.N, p.N]  = 0.0

p_autarky = Parameters(p.β, p.r, p.θ, p.δ, p.ρ, p.μ_z, p.ξ_z, p.ξ_ω, p.νᵈ, p.νᶠ, p.A, p.fᵈ, fᶠ_autarky, p.N)

# %% Autarky steady state
# With no foreign-born labor the domestic-born perform every task, so λ = 0,
# Z = E[z(τ)^b]^(1/b) = exp(μ_z + b·ξ_z²/2) and L = Lᵈ. With constant returns
# and capital rented at a fixed rate, the domestic wage is then the same
# constant in every US state, whatever its labor supply. The foreign-born are
# all in Rest of World, which is absorbing for them.
b = p.ρ / (1 - p.ρ)
κ = (p.θ / (p.r + p.δ))^(p.θ / (1 - p.θ))
Wᵈ_us_autarky = (1 - p.θ) * p.A * κ * exp(p.μ_z + 0.5 * b * p.ξ_z^2)

Wᵈ_autarky = vcat(fill(Wᵈ_us_autarky, L), Wᵈ_ROW_eff)
Vᵈ_autarky, Πᵈ_autarky = solve_value_choiceprobs(Wᵈ_autarky, p_autarky.fᵈ, p.β, p.νᵈ)
Lᵈ_autarky = stationary_distribution(Πᵈ_autarky, total_Lᵈ)

Lᶠ_autarky = vcat(zeros(L), total_Lᶠ)
Vᶠ_ROW_autarky = log(Wᶠ_ROW_eff) / (1 - p.β)

# %% Gains from migration: baseline relative to autarky
# Consumption equivalent g solves ln(1+g)/(1-β) = V_baseline - V_autarky: the
# permanent proportional change in consumption, valued at a home-bias weight
# of one, that is worth the same as moving from autarky to the baseline.
consumption_equivalent(V_baseline, V_autarky) = exp((1 - p.β) * (V_baseline - V_autarky)) - 1

wage_gainᵈ = ss.Wᵈ[1:L] ./ Wᵈ_us_autarky .- 1
ceᵈ = consumption_equivalent.(ss.Vᵈ[1:L], Vᵈ_autarky[1:L])
ceᶠ_ROW = consumption_equivalent(ss.Vᶠ[p.N], Vᶠ_ROW_autarky)

# Foreign-born living in US state l in the baseline are in Rest of World under
# autarky, so their comparison is the baseline value in l against the autarky
# value in Rest of World. The wage comparison uses the measured Rest of World
# wage; the values use the home-bias-adjusted one.
wage_gainᶠ = ss.Wᶠ[1:L] ./ Wᶠ_ROW .- 1
ceᶠ = consumption_equivalent.(ss.Vᶠ[1:L], Vᶠ_ROW_autarky)

states = load_bilateral_costs().states
weightsᵈ = ss.Lᵈ[1:L] ./ sum(ss.Lᵈ[1:L])
ceᵈ_mean = sum(weightsᵈ .* ceᵈ)
wage_gainᵈ_mean = sum(weightsᵈ .* wage_gainᵈ)
weightsᶠ = ss.Lᶠ[1:L] ./ sum(ss.Lᶠ[1:L])
ceᶠ_mean = sum(weightsᶠ .* ceᶠ)
wage_gainᶠ_mean = sum(weightsᶠ .* wage_gainᶠ)

println("Autarky domestic wage (all US states) = ", round(Wᵈ_us_autarky, digits = 1))
println("Domestic wage, baseline relative to autarky: mean ", round(100wage_gainᵈ_mean, digits = 3),
        "%  min ", round(100minimum(wage_gainᵈ), digits = 3), "%  max ", round(100maximum(wage_gainᵈ), digits = 3), "%")
println("Domestic consumption equivalent, baseline relative to autarky: mean ", round(100ceᵈ_mean, digits = 3),
        "%  min ", round(100minimum(ceᵈ), digits = 3), "%  max ", round(100maximum(ceᵈ), digits = 3), "%")
println("Foreign-born in Rest of World, consumption equivalent: ", round(100ceᶠ_ROW, digits = 5), "%")
println("Foreign-born in the US, wage relative to Rest of World: mean ", round(100wage_gainᶠ_mean, digits = 1),
        "%  min ", round(100minimum(wage_gainᶠ), digits = 1), "%  max ", round(100maximum(wage_gainᶠ), digits = 1), "%")
println("Foreign-born in the US, consumption equivalent: mean ", round(100ceᶠ_mean, digits = 2),
        "%  min ", round(100minimum(ceᶠ), digits = 2), "%  max ", round(100maximum(ceᶠ), digits = 2), "%")
println("Rest of World foreign wage: measured ", round(Wᶠ_ROW, digits = 1), "  home-bias adjusted ", round(Wᶠ_ROW_eff, digits = 1),
        "   baseline US foreign wage: mean ", round(sum(weightsᶠ .* ss.Wᶠ[1:L]), digits = 1))
println("Domestic-born share abroad: baseline ", round(100ss.Lᵈ[p.N] / total_Lᵈ, digits = 3),
        "%  autarky ", round(100Lᵈ_autarky[p.N] / total_Lᵈ, digits = 3), "%")

# %% Table: consumption-equivalent gains from migration, baseline relative to autarky (percent)
# The domestic-born mean weights states by baseline domestic-born labor.
using Printf

pct2(x) = @sprintf("%.2f", 100x)
open(joinpath(tables, "AutarkyGains.tex"), "w") do io
    println(io, "{")
    println(io, "\\begin{tabular}{lccc}")
    println(io, "\\toprule")
    println(io, "Consumption equivalent & Mean & Min & Max \\\\")
    println(io, "\\midrule")
    println(io, "Domestic-born, by US state & ", pct2(ceᵈ_mean), " & ", pct2(minimum(ceᵈ)), " & ", pct2(maximum(ceᵈ)), " \\\\")
    println(io, "Foreign-born in Rest of World & ", pct2(ceᶠ_ROW), " & & \\\\")
    println(io, "\\bottomrule")
    println(io, "\\end{tabular}")
    println(io, "}")
end
println("Wrote ", joinpath(tables, "AutarkyGains.tex"))

# %% Map: domestic-born consumption-equivalent gain by state (50 states, DC not drawn)
using CairoMakie, Shapefile

shp = Shapefile.Table(joinpath(data, "cb_2023_us_state_500k", "cb_2023_us_state_500k.shp"))
ce_by_fips = Dict(states[l] => 100 * ceᵈ[l] for l in 1:L if states[l] != "11")

"""Albers equal-area conic projection on the unit sphere (degrees in)."""
function albers(lon, lat; lat1, lat2, lat0, lon0)
    φ, λ, φ1, φ2, φ0, λ0 = deg2rad.((lat, lon, lat1, lat2, lat0, lon0))
    n  = (sin(φ1) + sin(φ2)) / 2
    C  = cos(φ1)^2 + 2n * sin(φ1)
    r  = sqrt(C - 2n * sin(φ)) / n
    r0 = sqrt(C - 2n * sin(φ0)) / n
    return Point2f(r * sin(n * (λ - λ0)), r0 - r * cos(n * (λ - λ0)))
end

"""
Projected outlines of one state: one vector of points per ring of the
boundary file, thinned to at most `maxpts` points, with rings shorter than
`minpts` (small islands) dropped. `place` rescales and moves the outline, used
to set Alaska and Hawaii below the contiguous states.
"""
function state_rings(geom; proj, place = identity, minpts = 60, maxpts = 400)
    pts    = geom.points
    starts = vcat(geom.parts .+ 1, length(pts) + 1)
    rings  = Vector{Vector{Point2f}}()
    for k in 1:length(starts) - 1
        idx = starts[k]:starts[k + 1] - 1
        length(idx) < minpts && continue
        step = max(1, length(idx) ÷ maxpts)
        # longitudes east of the antimeridian (western Aleutians) wrap to the west
        push!(rings, [place(proj(pts[i].x > 0 ? pts[i].x - 360 : pts[i].x, pts[i].y)) for i in idx[1:step:end]])
    end
    return rings
end

proj_main = (lon, lat) -> albers(lon, lat; lat1 = 29.5, lat2 = 45.5, lat0 = 37.5, lon0 = -96.0)
proj_ak   = (lon, lat) -> albers(lon, lat; lat1 = 55.0, lat2 = 65.0, lat0 = 62.0, lon0 = -152.0)
proj_hi   = (lon, lat) -> albers(lon, lat; lat1 = 8.0,  lat2 = 18.0, lat0 = 20.5, lon0 = -157.5)
place_ak  = q -> Point2f(0.35f0 * q[1] - 0.30f0, 0.35f0 * q[2] - 0.20f0)
place_hi  = q -> Point2f(q[1] - 0.17f0, q[2] - 0.215f0)

try
    set_theme!(fontsize = 16,
        fonts = (; regular = "Nimbus Roman", bold = "Nimbus Roman", italic = "Nimbus Roman"))
catch
    set_theme!(fontsize = 16)
end

crange  = extrema(values(ce_by_fips))
fig_map = Figure(size = (1000, 600))
ax_map  = Axis(fig_map[1, 1], aspect = DataAspect())

for (geom, fips) in zip(Shapefile.shapes(shp), shp.STATEFP)
    haskey(ce_by_fips, fips) || continue
    rings = fips == "02" ? state_rings(geom; proj = proj_ak, place = place_ak) :
            fips == "15" ? state_rings(geom; proj = proj_hi, place = place_hi) :
                           state_rings(geom; proj = proj_main)
    for ring in rings
        poly!(ax_map, ring; color = ce_by_fips[fips], colormap = :viridis, colorrange = crange,
            strokecolor = :white, strokewidth = 0.5)
    end
end

hidedecorations!(ax_map)
hidespines!(ax_map)
Colorbar(fig_map[1, 2]; colormap = :viridis, limits = crange,
    label = "Consumption-equivalent gain, domestic-born (%)")

outfile_map = joinpath(graphs, "AutarkyGainsMap.pdf")
save(outfile_map, fig_map)
println("Saved figure to $outfile_map")

# %% Bar chart: domestic-born consumption-equivalent gain by state, sorted (50 states, DC not drawn)
abbrev = Dict(zip(shp.STATEFP, shp.STUSPS))
order  = sort(collect(keys(ce_by_fips)), by = f -> ce_by_fips[f])

fig_bar = Figure(size = (1000, 450))
ax_bar  = Axis(fig_bar[1, 1], ylabel = "Consumption-equivalent gain, domestic-born (%)",
    xticks = (1:length(order), [abbrev[f] for f in order]), xticklabelrotation = π / 2, xticklabelsize = 11,
    xgridvisible = false, ygridvisible = false)
barplot!(ax_bar, 1:length(order), [ce_by_fips[f] for f in order], color = colorant"#1d3557")
hlines!(ax_bar, [100 * ceᵈ_mean], color = colorant"#8d99ae", linestyle = :dash, linewidth = 1.5, label = "Mean")
axislegend(ax_bar, position = :lt, framevisible = false)
hidespines!(ax_bar, :t, :r)
xlims!(ax_bar, 0.3, length(order) + 0.7)
ylims!(ax_bar, 0, nothing)

outfile_bar = joinpath(graphs, "AutarkyGainsBar.pdf")
save(outfile_bar, fig_bar)
println("Saved figure to $outfile_bar")
