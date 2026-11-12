# Builds output/tables/ParameterEstimates.tex (paper) and output/tables/ParameterEstimatesSlide.tex
# (slides; same two panels, with a text source in place of the paper's equation references).
# Panel A ("Estimating Equations") reports
# the parameters estimated directly from data under the log-normal
# task-assignment parameterization:
#   - (μ_z, ξ_ω, ξ_z), fit by EstimateCp.do's producer choice-probability
#     probit (eq:agg-foreign-cp)
#   - ρ, fit by AggSupply_Estimate.jl's factor-share NLS (eq:foreign-labor-share)
#   - (νᴰ, νᶠ), fit by EstimateScaleBetaAcs.do (eq:scale-estimating)
# The old Pareto/DFS parameterization's (θ, γ, μ) are gone -- γ, μ have no
# meaning outside that specification, and θ is no longer a single pooled
# scalar (it's state-varying, from CalibrateTheta.do's capital-FOC calibration,
# and doesn't fit this table's format). Panel B ("Moment Matching") reports parameters
# calibrated jointly against a set of target moments rather than read off a
# single closed-form estimating equation:
#   - (f^D_51, w^D_51, f^F_51, w^F_51, ι), the symmetric Rest of World migration
#     costs, the calibrated Rest of World wages and the quota on foreign-born
#     entry, from RowCostsEstimate.jl's baseline quota calibration
#     (RowCostsQuota.jld2), jointly matched to two stock shares and three flow
#     rates. The Target column pairs each parameter with the moment it most
#     directly governs: ι with the immigration rate (the realized inflow equals
#     the quota when it binds), the costs with the flow rates and the Rest of
#     World wages with the stock shares. The admission probability ξ is not a
#     parameter: it is solved from quota clearing, so it is left out of the table
#     and written as the LaTeX macro \xiSS to output/tables/CalibrationNumbers.tex.

using JLD2, Printf

include("Globals.jl")
include("Estimation_Funcs.jl")
include("RowCostsEstimate_Functions.jl")

cp        = load_cp_estimate()
aggsupply = load_aggsupply_estimate()
scale     = load_scale_estimate()
rowcosts  = load_row_costs_quota("baseline")

fmt(x) = @sprintf("%.3f", x)
pct(x) = @sprintf("%.3f\\%%", 100x)
# Dollar amounts and head counts with thousands separators
function thousands(x)
    d = string(round(Int, x))
    groups = [d[max(1, k - 2):k] for k in length(d):-3:1]
    return join(reverse(groups), "{,}")
end

# Elasticity of substitution between tasks implied by ρ, reported as its own row
σ̂ = 1 / (1 - aggsupply.ρ)

# (symbol, value, estimating-equation label, slide source)
panelA = [
    ("\\mu_z",         fmt(cp.μ_z),         "eq:agg-foreign-cp",      "Task shares (probit)"),
    ("\\xi_z",         fmt(cp.ξ_z),         "eq:agg-foreign-cp",      "Task shares (probit)"),
    ("\\xi_\\omega",   fmt(cp.ξ_ω),         "eq:agg-foreign-cp",      "Task shares (probit)"),
    ("\\rho",          fmt(aggsupply.ρ),    "eq:foreign-labor-share",
        "Foreign-born factor share \$\\rightarrow\\hat\\sigma = " * @sprintf("%.2f", σ̂) * "\$"),
    ("\\sigma",        fmt(σ̂),              "",                       "\$\\sigma = 1/(1-\\rho)\$"),
    ("\\nu^D",         fmt(scale.νᵈ),       "eq:scale-estimating",    "Migration flows on wages and future flows"),
    ("\\nu^F",         fmt(scale.νᶠ),       "eq:scale-estimating",    "Migration flows on wages and future flows"),
]

# (symbol, value, target title, model, data) -- shares/rates in percent so a
# tiny rate (the foreign immigration rate, ≈6.9e-5) doesn't round to 0.000 at
# three decimal places.
panelB_matched = [
    ("f^D_{51}", fmt(rowcosts.fᵈ_ROW), "Domestic Return Rate",
        pct(rowcosts.model_D_return), pct(rowcosts.target_D_return)),
    ("w^D_{51}", thousands(rowcosts.wᵈ_ROW), "Domestic Share Abroad",
        pct(rowcosts.model_D_abroad), pct(rowcosts.target_D_abroad)),
    ("f^F_{51}", fmt(rowcosts.fᶠ_ROW), "Foreign Application Share",
        pct(rowcosts.model_F_application), pct(rowcosts.target_F_application)),
    ("w^F_{51}", thousands(rowcosts.wᶠ_ROW), "Foreign Labor Share",
        pct(rowcosts.model_F_us), pct(rowcosts.target_F_us)),
    ("\\iota", thousands(rowcosts.ι_inflow_people), "Foreign Immigration Rate",
        pct(rowcosts.model_F_immigration), pct(rowcosts.target_F_immigration)),
]

"""
Write the two-panel parameter table. `source_header` is Panel A's column
heading and `source(row)` the text placed under it for each Panel A row.
"""
function write_parameter_table(outpath; source_header, source)

    io = IOBuffer()
    println(io, "{")
    println(io, "\\begin{tabular}{lclcc}")
    println(io, "\\toprule")
    println(io, " & Value & \\multicolumn{3}{l}{", source_header, "} \\\\")
    println(io, "\\midrule")
    println(io, "\\multicolumn{5}{l}{\\textit{Panel A. Estimating Equations}} \\\\[6pt]")
    for (i, row) in enumerate(panelA)
        sep = i == length(panelA) ? "\\\\[6pt]" : "\\\\"
        println(io, "\$", row[1], "\$ & ", row[2], " & \\multicolumn{3}{l}{", source(row), "} ", sep)
    end
    println(io, " & Value & Target & Model & Data \\\\")
    println(io, "\\midrule")
    println(io, "\\multicolumn{5}{l}{\\textit{Panel B. Moment Matching}} \\\\[6pt]")
    for (sym, val, target, model, dat) in panelB_matched
        println(io, "\$", sym, "\$ & ", val, " & ", target, " & ", model, " & ", dat, " \\\\")
    end
    println(io, "\\bottomrule")
    println(io, "\\end{tabular}")
    println(io, "}")

    write(outpath, String(take!(io)))
    println("Wrote ", outpath)

end

mkpath(tables)
write_parameter_table(joinpath(tables, "ParameterEstimates.tex");
    source_header = "Equation", source = row -> isempty(row[3]) ? row[4] : "Eq. \\ref{" * row[3] * "}")
write_parameter_table(joinpath(tables, "ParameterEstimatesSlide.tex");
    source_header = "Source", source = row -> row[4])

# Steady-state admission probability implied by quota clearing, for the text.
# Four decimals: at two it would print as 0.00.
open(joinpath(tables, "CalibrationNumbers.tex"), "w") do io
    println(io, "\\newcommand{\\xiSS}{", @sprintf("%.4f", rowcosts.ξ), "}")
    # Elasticity of substitution between tasks implied by the estimated ρ
    println(io, "\\newcommand{\\sigmaHat}{", @sprintf("%.2f", σ̂), "}")
end
println("Wrote ", joinpath(tables, "CalibrationNumbers.tex"))

# Definitions and data sources of the Panel B target moments (all 2015), for
# the paper's calibration discussion
moments = [
    ("Domestic share abroad",
     "US-born people living outside the US as a share of all US-born: those abroad plus US-born workers in the 50 states",
     "UN International Migrant Stock; ACS"),
    ("Domestic return rate",
     "US-born full-time workers aged 16 and older who lived abroad a year earlier, as a share of US-born people living abroad",
     "ACS 2016; UN International Migrant Stock"),
    ("Foreign labor share",
     "Foreign-born share of full-time workers in the 50 states",
     "ACS"),
    ("Foreign immigration rate",
     "Foreign-born full-time workers aged 16 and older who lived abroad a year earlier, as a share of the foreign-born population in Rest of World (world population less the US-born and the foreign-born in the US)",
     "ACS 2016; UN World Population Prospects; UN International Migrant Stock"),
    ("Foreign application share",
     "Adults naming the US as the country they would like to move to permanently, as a share of the same Rest of World population",
     "Gallup World Poll, 2015--2017 \\citep{esipova2018gallup}; UN World Population Prospects"),
]
open(joinpath(tables, "MomentDefinitions.tex"), "w") do io
    # Left-aligned paragraph columns: justified text spreads badly at these widths
    rag(w) = ">{\\raggedright\\arraybackslash}p{" * w * "\\linewidth}"
    println(io, "\\begin{tabular}{", rag("0.2"), rag("0.47"), rag("0.23"), "}")
    println(io, "\\toprule")
    println(io, "Moment & Description & Data Source \\\\")
    println(io, "\\midrule")
    for (name, desc, src) in moments
        println(io, name, " & ", desc, " & ", src, " \\\\[4pt]")
    end
    println(io, "\\bottomrule")
    println(io, "\\end{tabular}")
end
println("Wrote ", joinpath(tables, "MomentDefinitions.tex"))
