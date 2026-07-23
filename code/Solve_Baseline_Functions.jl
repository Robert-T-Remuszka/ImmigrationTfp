#================================================================
                        TYPE DECLARATIONS
================================================================#
"""
Parameters governing the baseline (factual, no-quota) sequential equilibrium.
Production parameters (ρ, θ, γ, μ) are the ones estimated in ProdFunc_Estimate.jl
"""
struct Parameters{T1 <: Real, T2 <: Integer}

    β::T1                                       # HH discount rate
    r::T1                                       # Capital rental rate
    δ::T1                                       # Capital depreciation rate
    ρ::T1                                       # CES parameter between foreign/domestic task aggregates
    θ::T1                                       # Capital share
    γ::T1                                       # Comparative-advantage schedule z(τ) = τ^γ
    μ::T1                                       # Pareto minimum of the variety distribution G(ω)
    ψ::T1                                       # AR(1) coefficient for the US-specific mobility cost mₜ
    νᵈ::T1                                      # Gumbel scale - domestic
    νᶠ::T1                                      # Gumbel scale - foreign

    N::T2                                       # Number of locations; convention is that region N is "Rest of World"
    wᵈ_row::T1                                  # Domestic wage earned in "Rest of World", real 1996 dollars per person
    wᶠ_row::T1                                  # Foreign wage earned in "Rest of World", real 1996 dollars per person
    Πᵈ₋::Matrix{T1}                             # Pre-period (1995→1996) choice probabilities, domestic (rows origin, cols destination)
    Πᶠ₋::Matrix{T1}                             # Pre-period choice probabilities, foreign
    Lᵈ₀::Vector{T1}                             # Time-zero (1996) domestic labor supplies, millions
    Lᶠ₀::Vector{T1}                             # Time-zero (1996) foreign labor supplies, millions
    Y₀::Vector{T1}                              # Time-zero (1996) real state GDP, millions of dollars (NaN for the Rest-of-World row)

end

"""
The baseline sequential equilibrium: wages, labor supplies, value changes and
choice probabilities over a T-period horizon. Πᵈ, Πᶠ hold [π₀, π₁, …, π_{T-2}] — length T-1 —
since a period-t choice probability governs the transition into period t+1.
"""
struct Soln{T1 <: Real, T2 <: Integer}

    Wᵈ::Matrix{T1}                              # Domestic wages (location × period)
    Wᶠ::Matrix{T1}                              # Foreign wages
    Lᵈ::Matrix{T1}                              # Domestic labor supplies
    Lᶠ::Matrix{T1}                              # Foreign labor supplies
    U̇ᵈ::Matrix{T1}                              # Domestic value changes
    U̇ᶠ::Matrix{T1}                              # Foreign value changes
    Πᵈ::Array{T1, 3}                            # Bilateral choice probabilities, domestic (origin × destination × period)
    Πᶠ::Array{T1, 3}                            # Bilateral choice probabilities, foreign

    T::T2                                       # Length of the solved horizon

end

#================================================================
                        DATA LOADING
================================================================#
"""
Load the pre-period (1995→1996) choice probabilities and time-zero (1996) labor stocks
from PiMat.dta. Rows and the pi_{F,D}_* columns are both ordered by FIPS code, with the
"Rest of World" region last — this is what gives Parameters struct its N-th-region convention.
"""
load_init_data() = DataFrame(load(joinpath(data, "PiMat.dta")))

"""
Load the current production-function estimate (ρ, θ, γ, μ) from ProductionFunction.jld2,
as fit by ProdFunc_Estimate.jl
"""
load_prodfunc_estimate() = load(joinpath(@__DIR__, "ProductionFunction.jld2"), "p_star")

"""
Load real 1996 state GDP (millions of dollars) from StateAnalysisPreTfp.dta, in the same
row order as Init_Data (matched on FIPS code — Init_Data's Origin column). The Rest-of-World
row has no GDP counterpart and is filled with NaN; it is never used (see solve_initial_wages).
"""
function load_gdp_1996(Init_Data::DataFrame)
    sa96 = DataFrame(load(joinpath(data, "StateAnalysisPreTfp.dta")))
    sa96 = sa96[sa96.year .== 1996, :]
    gdp_by_fips = Dict(sa96.statefip .=> Float64.(sa96.GDP))
    return [get(gdp_by_fips, o, NaN) for o in Init_Data.Origin]
end

#================================================================
                        CONSTRUCTOR FUNCTIONS
================================================================#
function Parameters(;
    β::T          = 0.96,
    r::T          = 1 / β - 1,
    δ::T          = 0.10,
    prodfunc::NamedTuple = load_prodfunc_estimate(),
    ρ::T          = prodfunc.ρ,
    θ::T          = prodfunc.θ,
    γ::T          = prodfunc.γ,
    μ::T          = prodfunc.μ,
    ψ::T          = 0.50,
    νᵈ::T         = 4.5,
    νᶠ::T         = 4.5,
    wᵈ_row::T     = 60077 / 0.77, # IRS Statistics of Income, US citizens abroad, 1996, deflated to 2009 USD
    wᶠ_row::T     = 4085  / 0.77, # World Bank GNI per capita, 1996, deflated to 2009 USD
    Init_Data::DataFrame = load_init_data()
    ) where {T <: Real}

    N = nrow(Init_Data)

    foreign_cols  = filter(c -> startswith(c, "pi_F_"), names(Init_Data))
    domestic_cols = filter(c -> startswith(c, "pi_D_"), names(Init_Data))
    Πᶠ₋ = Matrix{T}(Init_Data[:, foreign_cols])
    Πᵈ₋ = Matrix{T}(Init_Data[:, domestic_cols])

    Lᵈ₀ = Vector{T}(Init_Data[:, :Domestic_1996] ./ 1e+6)
    Lᶠ₀ = Vector{T}(Init_Data[:, :Foreign_1996]  ./ 1e+6)
    Y₀  = Vector{T}(load_gdp_1996(Init_Data))

    return Parameters(β, r, δ, ρ, θ, γ, μ, ψ, νᵈ, νᶠ, N, wᵈ_row, wᶠ_row, Πᵈ₋, Πᶠ₋, Lᵈ₀, Lᶠ₀, Y₀)

end

#================================================================
                    TEMPORARY EQUILIBRIUM 
================================================================#
"""
Task productivity Z and foreign-born task share λ implied by the relative wage w = wᴰ/wᶠ. 
Wrapper around ProdFunc.jl's TaskAggregates_μ.
"""
TaskAggregates(w; p::Parameters) = TaskAggregates_μ(p.ρ, p.γ, p.μ, w)

"""
Residual of the relative-wage / task-allocation condition:
    (wᴰ/wᶠ)^(-1/(1-ρ)) = (Lᴰ/Lᶠ) * λ/(1-λ)
"""
RelativeWageResidual(w, λ, lᵈ, lᶠ; ρ) = w^(-1 / (1 - ρ)) - (lᵈ / lᶠ) * (λ / (1 - λ))

"""
Scalar relative-wage residual (could be vectorized, but keep it simple), used by solve_scalar_wage below.
Solving for the level w_{l,t+1} directly given known w_{l,t} is the same equation as solving
for the proportional change ẇ_{l,t+1} = w_{l,t+1}/w_{l,t}
"""
function ScalarWageResidual(u, lᴰ_new, lᶠ_new; p::Parameters)

    w = exp(u)
    (; λ) = TaskAggregates(w; p)

    return RelativeWageResidual(w, λ, lᴰ_new, lᶠ_new; p.ρ)

end

"""
Newton solve for ScalarWageResidual. SolveTempEq
calls this once per (location, period) on every outer iteration of SolveBaseline, so its
per-call overhead dominates total runtime.
"""
function solve_scalar_wage(u0, lᴰ_new, lᶠ_new; p::Parameters, abstol::Real = 1e-6, maxiters::Integer = 100)

    u = u0
    for _ in 1:maxiters
        r = ScalarWageResidual(u, lᴰ_new, lᶠ_new; p)
        abs(r) < abstol && return u
        dr = ForwardDiff.derivative(x -> ScalarWageResidual(x, lᴰ_new, lᶠ_new; p), u)
        u -= r / dr
    end
    error("Scalar wage Newton solve failed to converge within $maxiters iterations")

end

"""
Solve the temporary equilibrium sequentially forward in time for each location, given the
already-updated labor stocks. Period 1 (t = 0, 1996) is fixed by the GDP anchor in
solve_initial_wages; periods 2, …, T are recovered from period 1's real-
dollar level via the differenced static system. TaskAggregates/LaborAggregate at the "old" level are cached from the previous
period's "new" level rather than recomputed — period t's new level is period t+1's old level.
The Rest-of-World wage (location N) is fixed exogenously.
"""
function SolveTempEq(S::Soln; p::Parameters, verbose::Bool = false)

    (; N, ρ) = p
    (; Wᵈ, Wᶠ, Lᵈ, Lᶠ, T) = S
    Wᵈ_new, Wᶠ_new = copy(Wᵈ), copy(Wᶠ)

    for l in 1:N - 1

        w_old  = Wᵈ_new[l, 1] / Wᶠ_new[l, 1]
        ta_old = TaskAggregates(w_old; p)
        L_old  = LaborAggregate(ta_old.λ, ρ, Lᶠ[l, 1], Lᵈ[l, 1])

        for t in 1:T - 1

            w_new = exp(solve_scalar_wage(log(w_old), Lᵈ[l, t + 1], Lᶠ[l, t + 1]; p))
            ta_new = TaskAggregates(w_new; p)
            L_new  = LaborAggregate(ta_new.λ, ρ, Lᶠ[l, t + 1], Lᵈ[l, t + 1])

            ratio = (ta_new.Z * L_new / (w_new * Lᵈ[l, t + 1] + Lᶠ[l, t + 1])) /
                    (ta_old.Z * L_old / (w_old * Lᵈ[l, t] + Lᶠ[l, t]))
            Wᶠ_new[l, t + 1] = Wᶠ_new[l, t] * ratio
            Wᵈ_new[l, t + 1] = w_new * Wᶠ_new[l, t + 1]

            w_old, ta_old, L_old = w_new, ta_new, L_new

        end

    end

    verbose && println("\tTemporary equilibrium solved for all locations and periods.")

    return Soln(Wᵈ_new, Wᶠ_new, S.Lᵈ, S.Lᶠ, S.U̇ᵈ, S.U̇ᶠ, S.Πᵈ, S.Πᶠ, S.T)

end

#================================================================
                INITIAL WAGE  (t = 0, real GDP anchor)
================================================================#
"""
Use resource feasibility, (1-θ)Y = wᴰlᴰ + wᶠlᶠ, using *observed* GDP Y in place
of the theoretical (θ/(r+δ))^{θ/(1-θ)}ZL construct in ResourceFeasResidual.

I use actual 1996 state GDP (Y₀, from StateAnalysisPreTfp.dta) at t = 0
only. This makes wᵈ_row, wᶠ_row (real dollars, from IRS SOI / World Bank) directly comparable
to the model's endogenous state wages without any separate calibration step.
"""
ResourceFeasResidualData(wᵈ, wᶠ, Y, lᵈ, lᶠ; θ) = (1 - θ) * Y / (wᵈ * lᵈ + wᶠ * lᶠ) - 1

"""
Bundles the relative-wage condition (3.7) and the data-anchored resource-feasibility residual
into a 2-vector, given u = [log(wᴰ/wᶠ), log(wᶠ)]. Used only at t = 0 — see solve_initial_wages.
"""
function InitialWageResidual(u, lᵈ, lᶠ, Y; p::Parameters)

    (; ρ, θ) = p
    w, wᶠ = exp(u[1]), exp(u[2])
    wᵈ    = w * wᶠ

    (; λ) = TaskAggregates(w; p)

    return [RelativeWageResidual(w, λ, lᵈ, lᶠ; ρ),
            ResourceFeasResidualData(wᵈ, wᶠ, Y, lᵈ, lᶠ; θ)]

end

"""
Backout init wages at t = 0 (1996) from actual state GDP (Y₀) and labor stocks
(Lᵈ₀, Lᶠ₀), via (3.7) + the data-anchored resource-feasibility residual above. The Rest-of-World
wage is exogenous and just carried over unchanged. This is the only place actual GDP data enters
the model; from t ≥ 1 the (forthcoming) proportional-change recursion takes over so that wage changes are
fully model-implied
"""
function solve_initial_wages(p::Parameters)

    (; N, Y₀, Lᵈ₀, Lᶠ₀, wᵈ_row, wᶠ_row) = p
    Wᵈ₀, Wᶠ₀ = zeros(N), zeros(N)
    Wᵈ₀[N], Wᶠ₀[N] = wᵈ_row, wᶠ_row

    for l in 1:N - 1

        u0 = [0.0, 0.0]
        f(u, _) = InitialWageResidual(u, Lᵈ₀[l], Lᶠ₀[l], Y₀[l]; p)
        sol = solve(NonlinearProblem(f, u0, nothing), NewtonRaphson(); maxiters = Int(1e6), abstol = 1e-6, reltol = 1e-6)
        sol.retcode == ReturnCode.Success || error("Initial wage bootstrap failed at location $l (retcode=$(sol.retcode))")
        Wᵈ₀[l] = exp(sol.u[1]) * exp(sol.u[2])
        Wᶠ₀[l] = exp(sol.u[2])

    end

    return Wᵈ₀, Wᶠ₀

end

#================================================================
                        SOLVE BASELINE ECONOMY
================================================================#
"""
N×N matrix of period-t mobility-cost multipliers. The US-specific mobility cost mₜ only
enters the migration problem (3.1) for moves from "Rest of World" (row N) into a US location
(columns 1:N-1); all other bilateral costs are unaffected.
"""
function cost_matrix(Ṁ_t, N)
    C = ones(N, N)
    C[N, 1:N - 1] .= Ṁ_t
    return C
end

"""
Update choice probabilities. Vectorized over (l, l') for each t: the
numerator matrix is the lagged probability times the forward value change and cost-change
terms, elementwise; rows are then normalized to sum to one.
"""
function UpdateChoiceProbabilities(S::Soln, Ṁ::Vector; p::Parameters)

    (; β, νᵈ, νᶠ, N, Πᵈ₋, Πᶠ₋) = p
    (; U̇ᵈ, U̇ᶠ, Πᵈ, Πᶠ, T) = S

    Πᵈ_new, Πᶠ_new = copy(Πᵈ), copy(Πᶠ)

    for t in 1:T - 1

        Πᵈ_lag = t == 1 ? Πᵈ₋ : Πᵈ[:, :, t - 1]
        Πᶠ_lag = t == 1 ? Πᶠ₋ : Πᶠ[:, :, t - 1]
        C      = cost_matrix(Ṁ[t], N)

        Numᵈ = Πᵈ_lag .* (U̇ᵈ[:, t + 1] .^ (β / νᵈ))' .* C .^ (-1 / νᵈ)
        Numᶠ = Πᶠ_lag .* (U̇ᶠ[:, t + 1] .^ (β / νᶠ))' .* C .^ (-1 / νᶠ)

        Πᵈ_new[:, :, t] = Numᵈ ./ sum(Numᵈ, dims = 2)
        Πᶠ_new[:, :, t] = Numᶠ ./ sum(Numᶠ, dims = 2)

    end

    return Soln(S.Wᵈ, S.Wᶠ, S.Lᵈ, S.Lᶠ, S.U̇ᵈ, S.U̇ᶠ, Πᵈ_new, Πᶠ_new, S.T)

end

"""
Update labor supplies, the factor-supply law of motion: Lₜ₊₁ = Πₜ' Lₜ.
"""
function UpdateLaborSupply(S::Soln; p::Parameters)

    (; N, Lᵈ₀, Lᶠ₀) = p
    (; Πᵈ, Πᶠ, T) = S

    Lᵈ_new, Lᶠ_new = zeros(N, T), zeros(N, T)
    Lᵈ_new[:, 1], Lᶠ_new[:, 1] = Lᵈ₀, Lᶠ₀

    for t in 1:T - 1
        Lᵈ_new[:, t + 1] = Πᵈ[:, :, t]' * Lᵈ_new[:, t]
        Lᶠ_new[:, t + 1] = Πᶠ[:, :, t]' * Lᶠ_new[:, t]
    end

    return Soln(S.Wᵈ, S.Wᶠ, Lᵈ_new, Lᶠ_new, S.U̇ᵈ, S.U̇ᶠ, S.Πᵈ, S.Πᶠ, S.T)

end

"""
Update value changes, by backward recursion from the boundary condition
U̇[:, T] = 1. Vectorized over l for each t: the inner sum over destinations l' is a
matrix-vector product of the (elementwise) lagged-probability/cost-change matrix against the
forward value changes.
"""
function UpdateValueChanges(S::Soln, Ṁ::Vector; p::Parameters)

    (; β, νᵈ, νᶠ, N, Πᵈ₋, Πᶠ₋) = p
    (; Πᵈ, Πᶠ, Wᵈ, Wᶠ, T) = S

    U̇ᵈ_new, U̇ᶠ_new = ones(N, T), ones(N, T)

    for t in T - 1:-1:1

        Πᵈ_lag = t == 1 ? Πᵈ₋ : Πᵈ[:, :, t - 1]
        Πᶠ_lag = t == 1 ? Πᶠ₋ : Πᶠ[:, :, t - 1]
        C      = cost_matrix(Ṁ[t], N)

        innerᵈ = (Πᵈ_lag .* C .^ (-1 / νᵈ)) * (U̇ᵈ_new[:, t + 1] .^ (β / νᵈ))
        innerᶠ = (Πᶠ_lag .* C .^ (-1 / νᶠ)) * (U̇ᶠ_new[:, t + 1] .^ (β / νᶠ))

        ẇᵈ = t == 1 ? ones(N) : Wᵈ[:, t] ./ Wᵈ[:, t - 1]
        ẇᶠ = t == 1 ? ones(N) : Wᶠ[:, t] ./ Wᶠ[:, t - 1]

        U̇ᵈ_new[:, t] = ẇᵈ .* innerᵈ .^ νᵈ
        U̇ᶠ_new[:, t] = ẇᶠ .* innerᶠ .^ νᶠ

    end

    return Soln(S.Wᵈ, S.Wᶠ, S.Lᵈ, S.Lᶠ, U̇ᵈ_new, U̇ᶠ_new, S.Πᵈ, S.Πᶠ, S.T)

end

#================================================================
                        MOBILITY-COST PATH
================================================================#
"""
The baseline economy: the US-specific mobility cost mₜ stays at its long-run mean forever, so
the proportional change Ṁₜ is 1 in every period. This is the reference equilibrium against
which future counterfactual cost experiments will be measured.
"""
no_shock() = t -> 1.0

#================================================================
                        OUTER SOLVER
================================================================#
"""
Take a damped convex combination of the newly updated U̇'s and their previous values, using
adaptive step size α: shrink α when the residual grows, grow it (up to 0.9) when it shrinks.
Returns the damped solution, the pre-damping residual, and the updated α.
"""
function damped_update(S::Soln, U̇ᵈ_prev, U̇ᶠ_prev, α::Real, prev_err::Real)

    err   = max(maximum(abs.(U̇ᵈ_prev .- S.U̇ᵈ)), maximum(abs.(U̇ᶠ_prev .- S.U̇ᶠ)))
    α_new = err > prev_err ? max(α * 0.5, 0.01) : min(α * 1.01, 0.9)

    U̇ᵈ_damped = α_new .* S.U̇ᵈ .+ (1 - α_new) .* U̇ᵈ_prev
    U̇ᶠ_damped = α_new .* S.U̇ᶠ .+ (1 - α_new) .* U̇ᶠ_prev
    S_damped = Soln(S.Wᵈ, S.Wᶠ, S.Lᵈ, S.Lᶠ, U̇ᵈ_damped, U̇ᶠ_damped, S.Πᵈ, S.Πᶠ, S.T)

    return S_damped, err, α_new

end

report_iteration(iter, err, α) = println("****************** ITERATION $iter COMPLETE: Outer err = $err (α = $α) *************************")

"""
Extend a converged Soln from T_old to T_new periods where T_old < T_new, holding
wages and choice probabilities fixed at their terminal values. Used to warm-start the solve
at a longer horizon once the current horizon isn't long enough to reach steady state — the
main speedup relative to solving a long horizon from scratch.
"""
function ExtendSoln(S::Soln, T_new::Integer)

    N, T_old = size(S.Wᵈ)
    @assert T_new > T_old "ExtendSoln requires T_new > T_old"
    dT = T_new - T_old

    Wᵈ_ext = hcat(S.Wᵈ, repeat(S.Wᵈ[:, end:end], 1, dT))
    Wᶠ_ext = hcat(S.Wᶠ, repeat(S.Wᶠ[:, end:end], 1, dT))
    U̇ᵈ_ext = hcat(S.U̇ᵈ, ones(N, dT))
    U̇ᶠ_ext = hcat(S.U̇ᶠ, ones(N, dT))
    Πᵈ_ext = cat(S.Πᵈ, repeat(S.Πᵈ[:, :, end:end], 1, 1, dT); dims = 3)
    Πᶠ_ext = cat(S.Πᶠ, repeat(S.Πᶠ[:, :, end:end], 1, 1, dT); dims = 3)

    Lᵈ_ext, Lᶠ_ext = hcat(S.Lᵈ, zeros(N, dT)), hcat(S.Lᶠ, zeros(N, dT))
    Πᵈ_term, Πᶠ_term = S.Πᵈ[:, :, end], S.Πᶠ[:, :, end]
    for t in T_old:T_new - 1
        Lᵈ_ext[:, t + 1] = Πᵈ_term' * Lᵈ_ext[:, t]
        Lᶠ_ext[:, t + 1] = Πᶠ_term' * Lᶠ_ext[:, t]
    end

    return Soln(Wᵈ_ext, Wᶠ_ext, Lᵈ_ext, Lᶠ_ext, U̇ᵈ_ext, U̇ᶠ_ext, Πᵈ_ext, Πᶠ_ext, T_new)

end

"""
Construct an initial guess for the baseline solution. Wages at t = 0 are the real-dollar GDP-
anchored from solve_initial_wages (fixed except for the exogenous Rest-of-World row);
all later periods start from that same level, choice probabilities are held at their pre-period
values, labor supplies are propagated forward under those probabilities, and value changes start
at their boundary value of 1.
"""
function Soln(p::Parameters; T0::Integer = 25)

    (; N, Lᵈ₀, Lᶠ₀, Πᵈ₋, Πᶠ₋) = p

    Wᵈ₀, Wᶠ₀ = solve_initial_wages(p)
    Wᵈ, Wᶠ   = repeat(Wᵈ₀, 1, T0), repeat(Wᶠ₀, 1, T0)

    Πᵈ, Πᶠ = repeat(Πᵈ₋, 1, 1, T0 - 1), repeat(Πᶠ₋, 1, 1, T0 - 1)

    Lᵈ, Lᶠ = zeros(N, T0), zeros(N, T0)
    Lᵈ[:, 1], Lᶠ[:, 1] = Lᵈ₀, Lᶠ₀
    for t in 1:T0 - 1
        Lᵈ[:, t + 1] = Πᵈ[:, :, t]' * Lᵈ[:, t]
        Lᶠ[:, t + 1] = Πᶠ[:, :, t]' * Lᶠ[:, t]
    end

    U̇ᵈ, U̇ᶠ = ones(N, T0), ones(N, T0)

    return Soln(Wᵈ, Wᶠ, Lᵈ, Lᶠ, U̇ᵈ, U̇ᶠ, Πᵈ, Πᶠ, T0)

end

"""
Solve for the baseline sequential competitive equilibrium: iterate {UpdateChoiceProbabilities, UpdateLaborSupply, SolveTempEq,
UpdateValueChanges} to convergence in U̇, extending the horizon whenever the solution hasn't
settled down to its steady-state boundary condition by period T-1.
"""
function SolveBaseline(p::Parameters = Parameters();
    T0::Integer = 25, outer_tol::Real = 1e-4, outer_maxiter::Integer = 10_000,
    init::Union{Soln, Nothing} = nothing, T_step::Integer = 10, ss_tol::Real = 1e-4,
    verbose::Bool = true)

    S   = isnothing(init) ? Soln(p; T0) : init
    Ṁ      = [no_shock()(t) for t in 1:S.T - 1]
    ss_err = 1 + ss_tol

    while ss_err >= ss_tol

        outer_err, outer_iter, prev_err, α = 1 + outer_tol, 0, Inf, 0.5

        while outer_err >= outer_tol && outer_iter <= outer_maxiter

            U̇ᵈ_prev, U̇ᶠ_prev = copy(S.U̇ᵈ), copy(S.U̇ᶠ)

            S = UpdateChoiceProbabilities(S, Ṁ; p)
            S = UpdateLaborSupply(S; p)
            S = SolveTempEq(S; p, verbose = verbose && outer_iter % 50 == 0)
            S = UpdateValueChanges(S, Ṁ; p)

            S, outer_err, α = damped_update(S, U̇ᵈ_prev, U̇ᶠ_prev, α, prev_err)
            prev_err   = outer_err
            outer_iter += 1

            verbose && outer_iter % 50 == 0 && report_iteration(outer_iter, outer_err, α)

        end

        # The T-th column of U̇ is already 1 by construction; check convergence to steady state
        # in the penultimate period instead.
        ss_err = max(maximum(abs.(S.U̇ᵈ[:, end - 1] .- 1)), maximum(abs.(S.U̇ᶠ[:, end - 1] .- 1)))

        if ss_err >= ss_tol
            verbose && println("Extending T from $(S.T) to $(S.T + T_step) (ss_err = $ss_err)")
            T_old = S.T
            S  = ExtendSoln(S, T_old + T_step)
            append!(Ṁ, [no_shock()(t) for t in T_old:T_old + T_step - 1])
        end

    end

    return S

end

#================================================================
                        DIAGNOSTICS
================================================================#
"""
Recover the converged task-productivity series Z(l, t) implied by a solved Soln.
"""
ComputeZ(S::Soln; p::Parameters) = getproperty.(TaskAggregates_μ.(p.ρ, p.γ, p.μ, S.Wᵈ ./ S.Wᶠ), :Z)