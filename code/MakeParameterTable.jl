using JLD2, Printf

include(joinpath(@__DIR__, "Globals.jl"))

"""
Structural parameter estimates table: production-function parameters (ρ, θ, γ, μ) from
ProductionFunction.jld2, indirect-inference parameters (νᴰ, νᶠ, ψ, σ) from
IndInfEstimate.jld2. Parameters go across columns (not the usual rows-of-coefficients
convention used for the regression tables in output/tables) since every entry here is a
single scalar, not a coefficient with a standard error.
"""

fmt(x) = @sprintf("%.3f", x)

prodfunc = load(joinpath(@__DIR__, "ProductionFunction.jld2"), "p_star")
indinf   = load(joinpath(data, "IndInfEstimate.jld2"))

symbols = ["\$\\rho\$", "\$\\theta\$", "\$\\gamma\$", "\$\\mu\$",
           "\$\\nu^D\$", "\$\\nu^F\$", "\$\\psi\$", "\$\\sigma\$"]

values = [prodfunc.ρ, prodfunc.θ, prodfunc.γ, prodfunc.μ,
          indinf["νᴰ̂"], indinf["νᶠ̂"], indinf["ψ̂"], indinf["σ̂"]]

group_row  = "            &\\multicolumn{4}{c}{Production Function}&\\multicolumn{4}{c}{Indirect Inference}\\\\\\cmidrule(lr){2-5}\\cmidrule(lr){6-9}"
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