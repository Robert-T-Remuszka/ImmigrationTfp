# Builds output/tables/ParameterEstimates.tex. Panel A ("Directly Inferred")
# reports the parameters estimated directly from data -- the production
# function (ρ, θ, γ, μ), fit by ProdFunc_Estimate.jl, and the migration
# elasticities (νᴰ, νᶠ), fit by EstimateScaleBetaAcs.do -- each row's Source
# column pointing to the paper's estimating equation. Panel B ("Indirectly
# Inferred") reports (ψ, σ), which will be estimated by indirect inference
# against the section-4 LPIV impulse responses; left blank until that step
# (and the accompanying Jacobian/target-moment table) exists -- see
# Estimation_Funcs.jl's load_irf_estimates() for the current state of that
# handoff.

using JLD2, Printf

include("Globals.jl")
include("Estimation_Funcs.jl")

prodfunc = load(joinpath(@__DIR__, "ProductionFunction.jld2"), "p_star")
scale    = load_scale_estimate()

fmt(x) = @sprintf("%.3f", x)

# (symbol, value, estimating-equation label)
panelA = [
    ("\\rho",   fmt(prodfunc.ρ), "eq:prod-fcn-est"),
    ("\\theta", fmt(prodfunc.θ), "eq:prod-fcn-est"),
    ("\\gamma", fmt(prodfunc.γ), "eq:prod-fcn-est"),
    ("\\mu",    fmt(prodfunc.μ), "eq:prod-fcn-est"),
    ("\\nu^D",  fmt(scale.νᵈ),   "eq:scale-estimating"),
    ("\\nu^F",  fmt(scale.νᶠ),   "eq:scale-estimating"),
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
