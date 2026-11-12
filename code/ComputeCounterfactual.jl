# %% Setup
using DataFrames, StatFiles, JLD2

include("Globals.jl")
include("TransitionSolve_Functions.jl")

#=================================================================
Autarky counterfactual, both steady state to steady state and along the
transition. The baseline is the calibrated non-stochastic steady state with
the quota ι on foreign-born entry. The counterfactual sets the quota to zero,
so the admission probability ξ is zero: the foreign-born already in the US
still draw taste shocks and leave, and none can enter, so the US has no
foreign-born labor in the new steady state.
=================================================================#

# %% Baseline: calibrated parameters and steady state
fixed = load_fixed_parameters()
cal = load_row_costs_quota("baseline")
p = Parameters(fixed, cal.fᵈ_ROW, cal.fᵈ_ROW, cal.fᶠ_ROW, cal.fᶠ_ROW)
wᵈ_ROW, wᶠ_ROW = cal.wᵈ_ROW, cal.wᶠ_ROW

# Measured Rest of World wages, for the wage comparisons only
(; Wᵈ_ROW, Wᶠ_ROW) = load_row_income()

tp = load_total_population()
total_Lᵈ = 1.0
total_Lᶠ = tp.total_Lᶠ / tp.total_Lᵈ

t0 = time()
ss = solve_steady_state(p; Wᵈ_ROW = wᵈ_ROW, Wᶠ_ROW = wᶠ_ROW, total_Lᵈ, total_Lᶠ,
    W0ᵈ = load_wages().Wᵈ, W0ᶠ = load_wages().Wᶠ, ι = cal.ι, ξ0 = cal.ξ)
println("Baseline steady state solved in ", round(time() - t0, digits = 1), "s, ξ = ", ss.ξ)

L = p.N - 1
println("Baseline foreign-born share of US labor = ", sum(ss.Lᶠ[1:L]) / (sum(ss.Lᶠ[1:L]) + sum(ss.Lᵈ[1:L])))

# %% Autarky steady state
# With no foreign-born labor the domestic-born perform every task, so λ = 0,
# Z = E[z(τ)^b]^(1/b) = exp(μ_z + b·ξ_z²/2) and L = Lᵈ. With constant returns
# and capital rented at a fixed rate, the domestic wage is then the same
# constant in every US state, whatever its labor supply. The foreign-born are
# all in Rest of World, which is absorbing for them since ξ = 0.
b = p.ρ / (1 - p.ρ)
κ = (p.θ / (p.r + p.δ))^(p.θ / (1 - p.θ))
Wᵈ_us_autarky = (1 - p.θ) * p.A * κ * exp(p.μ_z + 0.5 * b * p.ξ_z^2)

Wᵈ_autarky = vcat(fill(Wᵈ_us_autarky, L), wᵈ_ROW)
Vᵈ_autarky, Πᵈ_autarky = solve_value_choiceprobs(Wᵈ_autarky, p.fᵈ, p.β, p.νᵈ)
Lᵈ_autarky = stationary_distribution(Πᵈ_autarky, total_Lᵈ)

Lᶠ_autarky = vcat(zeros(L), total_Lᶠ)
Vᶠ_ROW_autarky = log(wᶠ_ROW) / (1 - p.β)

# %% Gains from migration: baseline relative to autarky
# Consumption equivalent g solves ln(1+g)/(1-β) = V_baseline - V_autarky: the
# permanent proportional change in consumption that is worth the same as
# moving from autarky to the baseline.
consumption_equivalent(V_baseline, V_autarky) = exp((1 - p.β) * (V_baseline - V_autarky)) - 1

wage_gainᵈ = ss.Wᵈ[1:L] ./ Wᵈ_us_autarky .- 1
ceᵈ = consumption_equivalent.(ss.Vᵈ[1:L], Vᵈ_autarky[1:L])
ceᶠ_ROW = consumption_equivalent(ss.Vᶠ[p.N], Vᶠ_ROW_autarky)

# Foreign-born living in US state l in the baseline are in Rest of World under
# autarky, so their comparison is the baseline value in l against the autarky
# value in Rest of World. The wage comparison uses the measured Rest of World
# wage; the values use the calibrated one.
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
println("Rest of World foreign wage: measured ", round(Wᶠ_ROW, digits = 1), "  calibrated ", round(wᶠ_ROW, digits = 1),
        "   baseline US foreign wage: mean ", round(sum(weightsᶠ .* ss.Wᶠ[1:L]), digits = 1))
println("Domestic-born share abroad: baseline ", round(100ss.Lᵈ[p.N] / total_Lᵈ, digits = 3),
        "%  autarky ", round(100Lᵈ_autarky[p.N] / total_Lᵈ, digits = 3), "%")

# %% Transition: the quota is set to zero unexpectedly at t=0
# Labor at t=0 is predetermined at the baseline steady state and the
# foreign-born already in the US leave only through their own choices.
T_transition = 300
path = solve_transition(p, ss, repeat(p.fᶠ, 1, 1, T_transition + 1); ι_path = zeros(T_transition + 1))

# Consumption equivalent ς along the transition: ln(1+ς)/(1-β) = V̂₀ - V, the
# permanent proportional change in baseline consumption worth the same as
# closing the border. ς<0 is a loss from closing it.
ςᵈ = consumption_equivalent.(path.Vᵈ[:, 1], ss.Vᵈ)
ςᶠ = consumption_equivalent.(path.Vᶠ[:, 1], ss.Vᶠ)

# Common ς whose log is the population-weighted mean of the log ς's
log_mean(μ, ς) = exp(sum(μ .* log.(1 .+ ς))) - 1
ςᵈ_mean = log_mean(weightsᵈ, ςᵈ[1:L])
ςᶠ_mean = log_mean(weightsᶠ, ςᶠ[1:L])

# Wage and option-value annuities in the value decomposition, by location.
# The option value is -ν ln π_ll, the log stay probability, except for the
# foreign-born in Rest of World under the lottery: their value is
# ln w + βV' - ξ·ν·ln π̃, with π̃ the probability of not applying (rejected
# applicants also stay, so the realized stay probability is not the one that
# enters).
discounts = p.β .^ (0:T_transition - 1)
wage_term(Ŵ, W) = [(1 - p.β) * sum(discounts .* log.(Ŵ[l, 1:T_transition] ./ W[l])) for l in 1:p.N]
option_term(Π̂, Π, ν) = [-(1 - p.β) * ν * sum(discounts .* log.(Π̂[l, l, 1:T_transition] ./ Π[l, l])) for l in 1:p.N]
wage_termᵈ, option_termᵈ = wage_term(path.Wᵈ, ss.Wᵈ), option_term(path.Πᵈ, ss.Πᵈ, p.νᵈ)
wage_termᶠ, option_termᶠ = wage_term(path.Wᶠ, ss.Wᶠ), option_term(path.Πᶠ, ss.Πᶠ, p.νᶠ)

ξ_log_π̃(ξ, V_next) = ξ == 0 ? 0.0 : ξ * log1p(-leave_probability(V_next, p.fᶠ, p.β, p.νᶠ, p.N))
lottery_flow = [ξ_log_π̃(path.ξ[k], path.Vᶠ[:, k + 1]) - ξ_log_π̃(ss.ξ, ss.Vᶠ) for k in 1:T_transition]
option_termᶠ[p.N] = -(1 - p.β) * p.νᶠ * sum(discounts .* lottery_flow)

# The value difference equals the two annuities plus the terminal value term
terminal_termᵈ = (1 - p.β) * p.β^T_transition .* (path.Vᵈ[:, end] .- ss.Vᵈ)
terminal_termᶠ = (1 - p.β) * p.β^T_transition .* (path.Vᶠ[:, end] .- ss.Vᶠ)
gapᵈ = maximum(abs, wage_termᵈ .+ option_termᵈ .+ terminal_termᵈ .- log.(1 .+ ςᵈ))
gapᶠ = maximum(abs, wage_termᶠ .+ option_termᶠ .+ terminal_termᶠ .- log.(1 .+ ςᶠ))
println("Value difference vs. wage + option + terminal terms, max gap: domestic ", gapᵈ, "  foreign ", gapᶠ)

println("Domestic-born in US states, ς: mean ", round(100ςᵈ_mean, digits = 3),
        "%  min ", round(100minimum(ςᵈ[1:L]), digits = 3), "%  max ", round(100maximum(ςᵈ[1:L]), digits = 3),
        "%  (wage ", round(100sum(weightsᵈ .* wage_termᵈ[1:L]), digits = 3), "%, option value ", round(100sum(weightsᵈ .* option_termᵈ[1:L]), digits = 3), "%)")
println("Foreign-born in US states, ς: mean ", round(100ςᶠ_mean, digits = 3),
        "%  min ", round(100minimum(ςᶠ[1:L]), digits = 3), "%  max ", round(100maximum(ςᶠ[1:L]), digits = 3),
        "%  (wage ", round(100sum(weightsᶠ .* wage_termᶠ[1:L]), digits = 3), "%, option value ", round(100sum(weightsᶠ .* option_termᶠ[1:L]), digits = 3), "%)")
println("Domestic-born in Rest of World, ς: ", round(100ςᵈ[p.N], digits = 3), "%")
println("Foreign-born in Rest of World, ς: ", round(100ςᶠ[p.N], digits = 3), "%")
println("Foreign-born in the US at t=T relative to the baseline: ", round(sum(path.Lᶠ[1:L, end]) / sum(ss.Lᶠ[1:L]), digits = 4))

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

# %% Bar chart: transition consumption equivalent by location, split into wage
# and option-value terms, for domestic-born and foreign-born workers: the 50
# states sorted by the total, then Rest of World (ROW) last. For the
# foreign-born in Rest of World, the potential immigrants, the option value is
# the lottery term computed above. Units are 100 × log points, i.e. percent of
# consumption. The sign is that of the gain from migration, the baseline
# relative to closing the border (minus ς), matching the impulse-response
# consumption equivalents in MakeIrfConsEq.jl, which are gains from an inflow.
function decomposition_panel!(ax, wage, option, total)
    order = vcat(sortperm(total[1:L]), p.N)
    x = 1:p.N
    barplot!(ax, repeat(x, 2), vcat(100 .* wage[order], 100 .* option[order]),
        dodge = repeat([1, 2], inner = p.N), color = repeat([colorant"#457b9d", colorant"#a8dadc"], inner = p.N))
    scatter!(ax, x, 100 .* total[order], color = colorant"#1d3557", markersize = 7)
    hlines!(ax, [0.0], color = :black, linewidth = 1)
    ax.xticks = (x, [l == p.N ? "ROW" : abbrev[states[l]] for l in order])
    hidespines!(ax, :t, :r)
    xlims!(ax, 0.3, p.N + 0.7)
end

fig_dec = Figure(size = (1200, 750))
ax_d = Axis(fig_dec[1, 1], ylabel = "Domestic-born (%)", xticklabelrotation = π / 2, xticklabelsize = 10,
    xgridvisible = false, ygridvisible = false)
ax_f = Axis(fig_dec[2, 1], ylabel = "Foreign-born (%)", xticklabelrotation = π / 2, xticklabelsize = 10,
    xgridvisible = false, ygridvisible = false)
decomposition_panel!(ax_d, -wage_termᵈ, -option_termᵈ, -log.(1 .+ ςᵈ))
decomposition_panel!(ax_f, -wage_termᶠ, -option_termᶠ, -log.(1 .+ ςᶠ))
Legend(fig_dec[1:2, 2],
    [PolyElement(color = colorant"#457b9d"), PolyElement(color = colorant"#a8dadc"),
     MarkerElement(marker = :circle, color = colorant"#1d3557", markersize = 9)],
    ["Wage", "Option value", "Total"], framevisible = false)

outfile_dec = joinpath(graphs, "AutarkyTransitionDecomp.pdf")
save(outfile_dec, fig_dec)
println("Saved figure to $outfile_dec")

# %% Headline numbers of the transition as LaTeX macros, for \input in the paper.
# Gains from migration in percent: the baseline relative to closing the border.
# The total is the consumption equivalent itself, 1/(1+ς) - 1, aggregated in
# levels; the wage and option-value terms are 100 × log points, as in the
# decomposition chart. Means over US states weight by baseline steady-state
# employment of the same nativity; Rest of World is reported on its own.
gainᵈ, gainᶠ = 1 ./ (1 .+ ςᵈ) .- 1, 1 ./ (1 .+ ςᶠ) .- 1
fmt(x) = (v = round(100x, digits = 2); @sprintf("%.2f", v == 0 ? 0.0 : v))
# Each number is also written as name*Abs, its absolute value, for text that
# states the direction in words ("lose").
texmacro(io, name, x) = (println(io, "\\newcommand{\\$(name)}{", fmt(x), "}");
                         println(io, "\\newcommand{\\$(name)Abs}{", fmt(abs(x)), "}"))
open(joinpath(tables, "AutarkyCeNumbers.tex"), "w") do io
    for (tag, gain, wage, option, w) in (("Dom", gainᵈ, -wage_termᵈ, -option_termᵈ, weightsᵈ),
                                         ("For", gainᶠ, -wage_termᶠ, -option_termᶠ, weightsᶠ))
        texmacro(io, "autCe$(tag)TotalMean", sum(w .* gain[1:L]))
        texmacro(io, "autCe$(tag)TotalMin", minimum(gain[1:L]))
        texmacro(io, "autCe$(tag)TotalMax", maximum(gain[1:L]))
        texmacro(io, "autCe$(tag)WageMean", sum(w .* wage[1:L]))
        texmacro(io, "autCe$(tag)OptionMean", sum(w .* option[1:L]))
        texmacro(io, "autCe$(tag)OptionMin", minimum(option[1:L]))
        texmacro(io, "autCe$(tag)OptionMax", maximum(option[1:L]))
        texmacro(io, "autCe$(tag)RowTotal", gain[p.N])
        texmacro(io, "autCe$(tag)RowWage", wage[p.N])
        texmacro(io, "autCe$(tag)RowOption", option[p.N])
    end
    # The same results stated from the closure's side, for the slides: closing
    # the border relative to the baseline: ς itself and the terms without the sign flip.
    for (tag, ς, wage, option, w) in (("Dom", ςᵈ, wage_termᵈ, option_termᵈ, weightsᵈ),
                                      ("For", ςᶠ, wage_termᶠ, option_termᶠ, weightsᶠ))
        texmacro(io, "bcCe$(tag)TotalMean", sum(w .* ς[1:L]))
        texmacro(io, "bcCe$(tag)WageMean", sum(w .* wage[1:L]))
        texmacro(io, "bcCe$(tag)OptionMean", sum(w .* option[1:L]))
        texmacro(io, "bcCe$(tag)RowTotal", ς[p.N])
    end
end
println("Wrote ", joinpath(tables, "AutarkyCeNumbers.tex"))

# %% Table: gains from migration along the transition (percent), the baseline
# relative to closing the border. The total is the consumption equivalent; the
# wage and option-value terms are 100 × log points, as in the decomposition
# chart. Incumbents in US states: means weighted by baseline steady-state
# employment of the same nativity, with the minimum and maximum of the total
# across states. Rest of World is a single location for each nativity.
# Numbers in math mode so negatives get a minus sign, not a hyphen
mfmt(x) = "\$" * fmt(x) * "\$"
function write_transition_table(path, totalᵈ, totalᶠ, wageᵈ, wageᶠ, optionᵈ, optionᶠ)
    open(path, "w") do io
        println(io, "{")
        # Full line width, with the extra space spread between the columns
        println(io, "\\begin{tabular*}{\\linewidth}{@{\\extracolsep{\\fill}}lccccc}")
        println(io, "\\toprule")
        println(io, " & \\multicolumn{3}{c}{Total} & \\multicolumn{2}{c}{Decomposition} \\\\")
        println(io, "\\cmidrule{2-4} \\cmidrule{5-6}")
        println(io, " & Mean & Min & Max & Wage & Option value \\\\")
        println(io, "\\midrule")
        for (lab, total, wage, option, w) in (("Domestic-born", totalᵈ, wageᵈ, optionᵈ, weightsᵈ),
                                              ("Foreign-born", totalᶠ, wageᶠ, optionᶠ, weightsᶠ))
            println(io, lab, " & ", mfmt(sum(w .* total[1:L])), " & ", mfmt(minimum(total[1:L])), " & ",
                    mfmt(maximum(total[1:L])), " & ", mfmt(sum(w .* wage[1:L])), " & ", mfmt(sum(w .* option[1:L])), " \\\\")
        end
        for (lab, total, option) in (("Domestic-born, Rest of World", totalᵈ, optionᵈ),
                                     ("Foreign-born, Rest of World", totalᶠ, optionᶠ))
            # The Rest of World wage is a calibrated constant, so its wage term is left blank
            println(io, lab, " & ", mfmt(total[p.N]), " & & & & ", mfmt(option[p.N]), " \\\\")
        end
        println(io, "\\bottomrule")
        println(io, "\\end{tabular*}")
        println(io, "}")
    end
    println("Wrote ", path)
end

# Paper: gains from migration, the baseline relative to closing the border
write_transition_table(joinpath(tables, "AutarkyTransitionGains.tex"),
                       gainᵈ, gainᶠ, -wage_termᵈ, -wage_termᶠ, -option_termᵈ, -option_termᶠ)
# Slides: the same results from the closure's side, closing the border
# relative to the baseline, so they read like the 1% immigration shock
write_transition_table(joinpath(tables, "BorderClosureGainsSlide.tex"),
                       ςᵈ, ςᶠ, wage_termᵈ, wage_termᶠ, option_termᵈ, option_termᶠ)
