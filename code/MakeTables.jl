# Builds output/tables/ParameterEstimates.tex. Panel A ("Directly Inferred")
# reports the parameters estimated directly from data under the log-normal
# task-assignment parameterization:
#   - (μ_z, ξ_ω, ξ_z), fit by EstimateCp.do's producer choice-probability
#     probit (eq:agg-foreign-cp)
#   - ρ, fit by AggSupply_Estimate.jl's factor-share NLS (eq:foreign-labor-share)
#   - (νᴰ, νᶠ), fit by EstimateScaleBetaAcs.do (eq:scale-estimating)
# The old Pareto/DFS parameterization's (θ, γ, μ) are gone -- γ, μ have no
# meaning outside that specification, and θ is no longer a single pooled
# scalar (it's state-varying, from CalibrateTheta.do's capital-FOC calibration,
# and doesn't fit this table's format). Panel B ("Indirectly Inferred") reports
# (ψ, σ), which will be estimated by indirect inference against the section-4
# LPIV impulse responses; left blank until that step (and the accompanying
# Jacobian/target-moment table) exists -- see Estimation_Funcs.jl's
# load_irf_estimates() for the current state of that handoff.

using JLD2, Printf

include("Globals.jl")
include("Estimation_Funcs.jl")

cp        = load_cp_estimate()
aggsupply = load_aggsupply_estimate()
scale     = load_scale_estimate()

fmt(x) = @sprintf("%.3f", x)

# (symbol, value, estimating-equation label)
panelA = [
    ("\\mu_z",         fmt(cp.μ_z),         "eq:agg-foreign-cp"),
    ("\\xi_z",         fmt(cp.ξ_z),         "eq:agg-foreign-cp"),
    ("\\xi_\\omega",   fmt(cp.ξ_ω),         "eq:agg-foreign-cp"),
    ("\\rho",          fmt(aggsupply.ρ),    "eq:foreign-labor-share"),
    ("\\nu^D",         fmt(scale.νᵈ),       "eq:scale-estimating"),
    ("\\nu^F",         fmt(scale.νᶠ),       "eq:scale-estimating"),
]

# (symbol,) -- Value/Target/Model/Data all pending indirect inference
panelB = ["\\psi", "\\sigma"]

io = IOBuffer()
println(io, "{")
println(io, "\\begin{tabular}{lcccc}")
println(io, "\\toprule")
println(io, " & Value & \\multicolumn{3}{l}{Estimating Equation} \\\\")
println(io, "\\midrule")
println(io, "\\multicolumn{5}{l}{\\textit{Panel A. Directly Inferred}} \\\\[6pt]")
for (i, (sym, val, eqlabel)) in enumerate(panelA)
    sep = i == length(panelA) ? "\\\\[6pt]" : "\\\\"
    println(io, "\$", sym, "\$ & ", val, " & \\multicolumn{3}{l}{Eq. \\ref{", eqlabel, "}} ", sep)
end
println(io, " & Value & Target & Model & Data \\\\")
println(io, "\\midrule")
println(io, "\\multicolumn{5}{l}{\\textit{Panel B. Indirectly Inferred}} \\\\[6pt]")
for sym in panelB
    println(io, "\$", sym, "\$ & -- & -- & -- & -- \\\\")
end
println(io, "\\bottomrule")
println(io, "\\end{tabular}")
println(io, "}")

mkpath(tables)
outpath = joinpath(tables, "ParameterEstimates.tex")
write(outpath, String(take!(io)))
println("Wrote ", outpath)
