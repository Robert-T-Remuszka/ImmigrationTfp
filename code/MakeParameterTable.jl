using JLD2, Printf

include(joinpath(@__DIR__, "Globals.jl"))

"""
Structural parameter estimates table: production-function parameters (ρ, θ, γ, μ) from
ProductionFunction.jld2, plus placeholder baseline values (νᴰ, νᶠ, ψ, σ) standing in for the
indirect-inference estimate while that estimation is still running. Parameters go across
columns (not the usual rows-of-coefficients convention used for the regression tables in
output/tables) since every entry here is a single scalar, not a coefficient with a standard
error.
"""

fmt(x) = @sprintf("%.3f", x)

prodfunc = load(joinpath(@__DIR__, "ProductionFunction.jld2"), "p_star")

# Placeholder parameters (not the final estimate) — the best point found so far by the
# parallel global search of the still-running indirect-inference job (3865225), targeting
# only the Z and L IRFs with equal weighting; local refinement (xNES) is still ongoing.
νᴰ̂, νᶠ̂, ψ̂, σ̂ = 1.4724, 1.2582, -0.8704, 1.8112

symbols = ["\$\\rho\$", "\$\\theta\$", "\$\\gamma\$", "\$\\mu\$",
           "\$\\nu^D\$", "\$\\nu^F\$", "\$\\psi\$", "\$\\sigma\$"]

values = [prodfunc.ρ, prodfunc.θ, prodfunc.γ, prodfunc.μ, νᴰ̂, νᶠ̂, ψ̂, σ̂]

group_row  = "            &\\multicolumn{4}{c}{Production Function}&\\multicolumn{4}{c}{Baseline}\\\\\\cmidrule(lr){2-5}\\cmidrule(lr){6-9}"
header_row = "            &" * join(symbols, "&") * "\\\\"
value_row  = "            &" * join(fmt.(values), "&") * "\\\\"

table = """
{
\\begin{tabular}{l*{8}{c}}
\\toprule
$group_row
$header_row
\\midrule
$value_row
\\bottomrule
\\end{tabular}
}
"""

outfile = joinpath(tables, "ParameterEstimates.tex")
write(outfile, table)
println("Wrote parameter table to $outfile")
println(table)