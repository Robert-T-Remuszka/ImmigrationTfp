# %% Setup
using CSV, DataFrames, StatFiles, Printf, CairoMakie

include("Globals.jl")
include("SsSolve_Functions.jl")

#=================================================================
Consumption equivalents of a migration shock for incumbent workers, from the
local-projection IV impulse responses (MakeIrf.do's percent-change
specification, exported to LpivIrfPctChg_1990.csv). The first-order
approximation of the consumption equivalent, truncated at the last estimated
horizon H, is

    ς^n ≈ (1-β) Σ_{h=0}^{H} β^h [ IRF of ln w^n at h  -  ν^n · IRF of ln π^n_ll at h ],

for nativity n, with the wage and stay-probability responses of that
nativity. Each impulse response is per unit change in the log foreign-born
stock on impact, so ς^n is the percent change in consumption per one percent
increase in the foreign-born stock. Horizons 0..H carry 1-β^(H+1) of the
annuity weight, so these are the gains realized over that window.
=================================================================#

# %% Inputs
irf = CSV.read(joinpath(data, "LpivIrfPctChg_1990.csv"), DataFrame)
β = load_fixed_parameters().β
(; νᵈ, νᶠ) = load_scale_estimate()
H = maximum(irf.h)
discounts = β .^ (0:H)

irf_path(v) = sort(irf[irf.variable .== v, :], :h).beta

# %% Wage and option-value terms of the linear consumption equivalent
function linear_ce(wage, stay, ν)
    wage_part   = (1 - β) * sum(discounts .* irf_path(wage))
    option_part = -(1 - β) * ν * sum(discounts .* irf_path(stay))
    return (; wage_part, option_part, total = wage_part + option_part)
end

ceᵈ = linear_ce("Wage_Domestic_Resid", "Stay_Domestic", νᵈ)
ceᶠ = linear_ce("Wage_Foreign_Resid", "Stay_Foreign", νᶠ)

println("β = $β  νᵈ = $νᵈ  νᶠ = $νᶠ  horizons 0..$H carry ", round(100 * (1 - β^(H + 1)), digits = 1), "% of the annuity weight")
for (lab, c) in (("Domestic-born", ceᵈ), ("Foreign-born", ceᶠ))
    @printf("%-14s  wage %7.4f  option value %7.4f  total %7.4f  (percent per 1%% increase in the foreign-born stock)\n",
            lab, c.wage_part, c.option_part, c.total)
end

# %% Headline numbers as LaTeX macros, percent per one percent increase in the
# foreign-born stock, for \input in the paper, rounded to two decimals
# Each number is also written as name*Abs, its absolute value, for text that
# states the direction in words ("falls by").
fmt(x) = (v = round(x, digits = 2); @sprintf("%.2f", v == 0 ? 0.0 : v))
texmacro(io, name, x) = (println(io, "\\newcommand{\\$(name)}{", fmt(x), "}");
                         println(io, "\\newcommand{\\$(name)Abs}{", fmt(abs(x)), "}"))
open(joinpath(tables, "IrfCeNumbers.tex"), "w") do io
    for (tag, c) in (("Dom", ceᵈ), ("For", ceᶠ))
        texmacro(io, "irfCe$(tag)Total", c.total)
        texmacro(io, "irfCe$(tag)Wage", c.wage_part)
        texmacro(io, "irfCe$(tag)Option", c.option_part)
    end
    # Impulse responses at the last horizon (h = H, the tenth year), percent per
    # one percent increase in the foreign-born stock, with standard errors
    names = Dict("L" => "LaborAgg", "Supply_Foreign" => "LaborFor", "Supply_Domestic" => "LaborDom",
                 "Z" => "Productivity", "Stay_Domestic" => "StayDom", "Stay_Foreign" => "StayFor",
                 "Wage_Domestic_Resid" => "WageDom", "Wage_Foreign_Resid" => "WageFor", "Wage_Relative" => "WageRel")
    for r in eachrow(irf[irf.h .== H, :])
        texmacro(io, "irfTenYr$(names[r.variable])", r.beta)
        println(io, "\\newcommand{\\irfTenYr$(names[r.variable])Se}{", fmt(r.se), "}")
    end
end
println("Wrote ", joinpath(tables, "IrfCeNumbers.tex"))

# %% Bar chart: wage and option-value parts and the total, by nativity
try
    set_theme!(fontsize = 16,
        fonts = (; regular = "Nimbus Roman", bold = "Nimbus Roman", italic = "Nimbus Roman"))
catch
    set_theme!(fontsize = 16)
end

parts   = ["Wage", "Option value", "Total"]
colors  = [colorant"#457b9d", colorant"#a8dadc", colorant"#1d3557"]
values  = [ceᵈ.wage_part, ceᵈ.option_part, ceᵈ.total, ceᶠ.wage_part, ceᶠ.option_part, ceᶠ.total]
groups  = repeat([1, 2], inner = 3)
dodges  = repeat(1:3, outer = 2)

fig = Figure(size = (700, 450))
ax  = Axis(fig[1, 1], ylabel = "Consumption equivalent (%)",
    xticks = (1:2, ["Domestic-born", "Foreign-born"]), xgridvisible = false, ygridvisible = false)
barplot!(ax, groups, values, dodge = dodges, color = colors[dodges])
hlines!(ax, [0.0], color = :black, linewidth = 1)
Legend(fig[1, 2], [PolyElement(color = c) for c in colors], parts, framevisible = false)
hidespines!(ax, :t, :r)

outfile = joinpath(graphs, "IrfConsumptionEquivalents.pdf")
save(outfile, fig)
println("Saved figure to $outfile")
