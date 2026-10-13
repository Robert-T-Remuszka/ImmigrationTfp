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
#   - (f^D_ROW, B^D_ROW, f^F_ROW, B^F_ROW), the symmetric ROW migration costs
#     and nativity-specific home-bias amenities from RowCostsEstimate.jl,
#     jointly matched to two stock shares and two flow rates (see that
#     script's docstrings for why a symmetric cost alone can't do this
#     without the amenity). The Target column's parameter-to-moment pairing
#     is the one-to-one assignment maximizing total matched elasticity
#     (computed once via finite differences at the calibrated point; not
#     recomputed here since it doesn't change unless the calibration is
#     rerun) -- not just each parameter's individually-largest elasticity,
#     which for the domestic pair would double-assign both parameters to the
#     same moment.

using JLD2, Printf

include("Globals.jl")
include("Estimation_Funcs.jl")
include("RowCostsEstimate_Functions.jl")

cp        = load_cp_estimate()
aggsupply = load_aggsupply_estimate()
scale     = load_scale_estimate()
rowcosts  = load_row_costs()

fmt(x) = @sprintf("%.3f", x)
pct(x) = @sprintf("%.3f\\%%", 100x)

# Elasticity of substitution implied by ρ, shown next to ρ's source on the slide table
σ̂ = 1 / (1 - aggsupply.ρ)

# (symbol, value, estimating-equation label, slide source)
panelA = [
    ("\\mu_z",         fmt(cp.μ_z),         "eq:agg-foreign-cp",      "Task shares (probit)"),
    ("\\xi_z",         fmt(cp.ξ_z),         "eq:agg-foreign-cp",      "Task shares (probit)"),
    ("\\xi_\\omega",   fmt(cp.ξ_ω),         "eq:agg-foreign-cp",      "Task shares (probit)"),
    ("\\rho",          fmt(aggsupply.ρ),    "eq:foreign-labor-share",
        "Foreign-born factor share \$\\rightarrow\\hat\\sigma = " * @sprintf("%.2f", σ̂) * "\$"),
    ("\\nu^D",         fmt(scale.νᵈ),       "eq:scale-estimating",    "Migration flows on wages and future flows"),
    ("\\nu^F",         fmt(scale.νᶠ),       "eq:scale-estimating",    "Migration flows on wages and future flows"),
]

# (symbol, value, target title, model, data) -- shares/rates in percent so a
# tiny rate (the foreign immigration rate, ≈6.9e-5) doesn't round to 0.000 at
# three decimal places.
panelB_matched = [
    ("f^D_{\\text{ROW}}", fmt(rowcosts.fᵈ_ROW), "Domestic Return Rate",
        pct(rowcosts.model_D_return), pct(rowcosts.target_D_return)),
    ("B^D_{\\text{ROW}}", fmt(rowcosts.Bᵈ_ROW), "Domestic Share Abroad",
        pct(rowcosts.model_D_abroad), pct(rowcosts.target_D_abroad)),
    ("f^F_{\\text{ROW}}", fmt(rowcosts.fᶠ_ROW), "Foreign Immigration Rate",
        pct(rowcosts.model_F_immigration), pct(rowcosts.target_F_immigration)),
    ("B^F_{\\text{ROW}}", fmt(rowcosts.Bᶠ_ROW), "Foreign Labor Share",
        pct(rowcosts.model_F_us), pct(rowcosts.target_F_us)),
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
    source_header = "Equation", source = row -> "Eq. \\ref{" * row[3] * "}")
write_parameter_table(joinpath(tables, "ParameterEstimatesSlide.tex");
    source_header = "Source", source = row -> row[4])
