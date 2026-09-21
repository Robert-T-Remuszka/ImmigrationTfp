using DataFrames, StatFiles, JLD2, LinearAlgebra, ForwardDiff, NonlinearSolve

include("Estimation_Funcs.jl")
include("AggSupply_Functions.jl")

# %% Structs =================================================================

"""
Model parameters. fᵈ, fᶠ are N×N bilateral migration cost matrices, N = 51 US
states + Rest of World (the last row/column). A is the output scale
parameter (paper's Technology section) anchoring the model to real 2015
dollars.
"""
struct Parameters{T <: Real}
    β::T
    r::T
    θ::T
    δ::T
    ρ::T
    μ_z::T
    ξ_z::T
    ξ_ω::T
    νᵈ::T
    νᶠ::T
    A::T
    fᵈ::Matrix{T}
    fᶠ::Matrix{T}
    N::Int
end

"""
Steady-state values, wages, labor supplies, and choice probabilities, by
nativity. Wᵈ, Wᶠ are inputs here, not yet solved for by market clearing.
"""
struct Solution{T <: Real}
    Vᵈ::Vector{T}
    Vᶠ::Vector{T}
    Wᵈ::Vector{T}
    Wᶠ::Vector{T}
    Lᵈ::Vector{T}
    Lᶠ::Vector{T}
    Πᵈ::Matrix{T}
    Πᶠ::Matrix{T}
end

# %% Solver Functions =========================================================

softmax(x) = (e = exp.(x .- maximum(x)); e ./ sum(e))

"""
EMAX operator fixed point for one nativity, found via Newton's method on
g(V) = T(V) - V. The EMAX operator's Jacobian is β·Π (Π the choice
probabilities at the current V), so each step solves the N×N linear system
(I - βΠ)ΔV = g(V) -- quadratic convergence, no need for the many sweeps plain
fixed-point iteration would take.
"""
function solve_value_choiceprobs(W::Vector{T}, f::Matrix{T}, β::T, ν::T;
                                  tol::Real = 1e-10, maxiter::Integer = 100) where {T <: Real}

    N = length(W)
    lnW = log.(W)
    V = copy(lnW)
    Tv = similar(V)
    Π = zeros(T, N, N)

    for _ in 1:maxiter
        for l in 1:N
            x = (β .* V .- f[l, :]) ./ ν
            Π[l, :] .= softmax(x)
            Tv[l] = lnW[l] + ν * logsumexp(x)
        end
        g = Tv .- V
        maximum(abs.(g)) < tol && break
        V .+= (I - β .* Π) \ g
    end

    return V, Π

end

"""
Probability that a resident of location l leaves it, i.e. 1-Π[l,l], computed
directly from the softmax's off-diagonal weights rather than as one minus the
stored diagonal entry. When leaving l is very costly, Π[l,l] rounds to
exactly 1.0 in floating point -- the true leaving probability is smaller than
machine epsilon relative to 1 -- so 1-Π[l,l] silently returns exactly 0
regardless of how large the cost actually is, masking any further increase in
cost as "no additional effect" on a moment defined this way. Summing the
off-diagonal softmax weights directly avoids ever subtracting two near-equal
numbers, so it stays accurate all the way down to the smallest doubles.
"""
function leave_probability(V::Vector{T}, f::Matrix{T}, β::T, ν::T, l::Integer) where {T <: Real}

    x = (β .* V .- f[l, :]) ./ ν
    e = exp.(x .- maximum(x))

    return (sum(e) - e[l]) / sum(e)

end

"""
Stationary distribution of a row-stochastic transition matrix Π, scaled to
sum to `total`, via the GTH (Grassmann-Taksar-Heyman) algorithm. GTH
eliminates states one at a time from N down to 1, at each step redistributing
the eliminated state's outgoing probability back onto the remaining states --
every intermediate quantity is a sum, product, or ratio of nonnegative
numbers, so it cannot produce a negative component the way a generic linear
solve on (Π'-I) can. That matters here specifically: large migration costs
make Π nearly block-diagonal (nearly reducible), which is exactly where a
generic solve loses accuracy to floating-point cancellation and can return
small negative entries -- entries that are never legitimate, since Π (a
softmax output) has full support, so by Perron-Frobenius the true stationary
distribution is strictly positive everywhere.
"""
function stationary_distribution(Π::Matrix{T}, total::Real) where {T <: Real}

    N = size(Π, 1)
    a = Matrix{Float64}(Π) - I  # generator Q = Π - I: rows sum to zero
    col_save = zeros(Float64, N, N)
    S = zeros(Float64, N)

    for k in N:-1:2
        S[k] = sum(@view a[k, 1:k - 1])
        for i in 1:k - 1
            col_save[i, k] = a[i, k]
        end
        for i in 1:k - 1, j in 1:k - 1
            i == j && continue
            a[i, j] += a[i, k] * a[k, j] / S[k]
        end
    end

    π = zeros(Float64, N)
    π[1] = 1.0
    for k in 2:N
        π[k] = sum(π[i] * col_save[i, k] for i in 1:k - 1) / S[k]
    end

    return T.(π ./ sum(π) .* total)

end

"""
Residual of the relative-wage/task-allocation condition: w^(-1/(1-ρ)) =
(Lᵈ/Lᶠ)·(λ/(1-λ)), where λ is the foreign task share TaskAggregates_LN
implies at relative wage w = wᵈ/wᶠ.
"""
relative_wage_residual(w, λ, Lᵈ, Lᶠ, ρ) = w^(-1 / (1 - ρ)) - (Lᵈ / Lᶠ) * (λ / (1 - λ))

"""
Newton solve (in log w) for the relative wage w = wᵈ/wᶠ clearing the
task-allocation margin in one location, given its labor stocks.
"""
function solve_relative_wage(Lᵈ::Real, Lᶠ::Real, p::Parameters;
                              u0::Real = 0.0, tol::Real = 1e-10, maxiter::Integer = 100)

    resid(u) = relative_wage_residual(exp(u), TaskAggregates_LN(p.ρ, p.μ_z, p.ξ_ω, p.ξ_z, exp(u)).λ, Lᵈ, Lᶠ, p.ρ)

    u = u0
    for _ in 1:maxiter
        r = resid(u)
        abs(r) < tol && return exp(u)
        u -= r / ForwardDiff.derivative(resid, u)
    end
    error("Relative wage Newton solve failed to converge")

end

"""
Market-clearing wages Wᵈ, Wᶠ in one US location given its labor stocks
Lᵈ, Lᶠ: solve the relative wage (above), then recover levels in closed form
from resource feasibility with capital substituted out via the capital FOC,
(1-θ)Y = wᵈLᵈ + wᶠLᶠ where Y = A·(θ/(r+δ))^(θ/(1-θ))·Z·L.
"""
function market_clearing_wage(Lᵈ::Real, Lᶠ::Real, p::Parameters)

    w = solve_relative_wage(Lᵈ, Lᶠ, p)
    (; Z, λ, oneMinusλ) = TaskAggregates_LN(p.ρ, p.μ_z, p.ξ_ω, p.ξ_z, w)
    L = LaborAggregate(λ, p.ρ, Lᶠ, Lᵈ, oneMinusλ)
    Y = p.A * (p.θ / (p.r + p.δ))^(p.θ / (1 - p.θ)) * Z * L

    Wᶠ = (1 - p.θ) * Y / (w * Lᵈ + Lᶠ)
    Wᵈ = w * Wᶠ

    return Wᵈ, Wᶠ

end

"""
Residual of the full steady-state wage system: given log-wages for the 51 US
states (both nativities; Rest of World's wages are fixed exogenously), solve
the migration block for the labor allocation those wages generate, then
return the gap between the guessed wages and the wages market clearing
implies at that allocation. This is zero exactly at the non-stochastic
steady state.
"""
function steady_state_residual(logW::Vector{T}, params) where {T <: Real}

    (; p, Wᵈ_ROW, Wᶠ_ROW, total_Lᵈ, total_Lᶠ) = params
    L = p.N - 1

    Wᵈ = vcat(exp.(logW[1:L]), Wᵈ_ROW)
    Wᶠ = vcat(exp.(logW[L + 1:2L]), Wᶠ_ROW)

    s = Solution(p, Wᵈ, Wᶠ; total_Lᵈ, total_Lᶠ)

    res = similar(logW)
    for l in 1:L
        Wᵈ_new, Wᶠ_new = market_clearing_wage(s.Lᵈ[l], s.Lᶠ[l], p)
        res[l]     = log(Wᵈ_new) - logW[l]
        res[L + l] = log(Wᶠ_new) - logW[L + l]
    end

    return res

end

"""
Solve the full non-stochastic steady state: wages, values, labor supplies,
and choice probabilities jointly consistent with the migration block and
market clearing. Uses a finite-difference Jacobian rather than
differentiating through the nested EMAX/relative-wage Newton solves inside
the residual -- more robust, and cheap enough at this problem size (102
wage unknowns) that there's no real efficiency case for the riskier route.
"""
function solve_steady_state(p::Parameters; Wᵈ_ROW::Real, Wᶠ_ROW::Real,
                             total_Lᵈ::Real, total_Lᶠ::Real,
                             W0ᵈ::Vector{<:Real}, W0ᶠ::Vector{<:Real},
                             abstol::Real = 1e-8, maxiters::Integer = 1000)

    params = (; p, Wᵈ_ROW = Float64(Wᵈ_ROW), Wᶠ_ROW = Float64(Wᶠ_ROW), total_Lᵈ, total_Lᶠ)
    logW0 = log.(Float64.(vcat(W0ᵈ, W0ᶠ)))

    prob = NonlinearProblem(steady_state_residual, logW0, params)
    sol  = solve(prob, NewtonRaphson(autodiff = AutoFiniteDiff()); abstol, maxiters)

    sol.retcode == ReturnCode.Success || error("Steady-state solve failed (retcode=$(sol.retcode))")

    L  = p.N - 1
    Wᵈ = vcat(exp.(sol.u[1:L]), Wᵈ_ROW)
    Wᶠ = vcat(exp.(sol.u[L + 1:2L]), Wᶠ_ROW)

    return Solution(p, Wᵈ, Wᶠ; total_Lᵈ, total_Lᶠ)

end

# %% Helper Functions =========================================================

"""
Load the 51×51 fᵈ, fᶠ bilateral cost matrices from EstimateBilateralCosts.do's
saved BilateralCosts2015.dta, ordered by sorted state FIPS code.
"""
function load_bilateral_costs()

    df = DataFrame(load(joinpath(data, "BilateralCosts2015.dta")))
    states = sort(unique(df.Origin))
    idx = Dict(s => i for (i, s) in enumerate(states))
    L = length(states)

    fᵈ = zeros(L, L)
    fᶠ = zeros(L, L)
    for row in eachrow(df)
        i, j = idx[row.Origin], idx[row.Destination]
        fᵈ[i, j] = row.f_Domestic
        fᶠ[i, j] = row.f_Foreign
    end

    return (; fᵈ, fᶠ, states)

end

"""
Load every Parameters input except the four ROW-linkage costs, once. Splitting
this disk-reading step out of Parameters() means a grid or threaded search
over those four costs can build thousands of Parameters objects without
touching disk on the hot path -- JLD2 reads aren't thread-safe for concurrent
access, and re-reading the same five files on every evaluation is wasted I/O
regardless of threading.
"""
function load_fixed_parameters(; β::Real = 0.96)

    r  = 1 / β - 1
    θδ = load_theta_delta()
    θ, δ = θδ.theta_l[1], θδ.delta_l[1]
    ρ = load_aggsupply_estimate().ρ
    (; μ_z, ξ_ω, ξ_z) = load_cp_estimate()
    (; νᵈ, νᶠ) = load_scale_estimate()
    (; fᵈ, fᶠ) = load_bilateral_costs()
    A = calibrate_scale(θ, δ, r, ρ, μ_z, ξ_ω, ξ_z)
    L = size(fᵈ, 1)

    return (β = Float64(β), r = Float64(r), θ = Float64(θ), δ = Float64(δ), ρ = Float64(ρ),
            μ_z = Float64(μ_z), ξ_z = Float64(ξ_z), ξ_ω = Float64(ξ_ω),
            νᵈ = Float64(νᵈ), νᶠ = Float64(νᶠ), A = Float64(A),
            fᵈ = Matrix{Float64}(fᵈ), fᶠ = Matrix{Float64}(fᶠ), L = L, N = L + 1)

end

"""
Build Parameters from pre-loaded `fixed` (see load_fixed_parameters) plus the
four ROW-linkage costs, touching no disk. The 51×51 interior keeps its
estimated f^n_{ll'}=f^n_{l'l} symmetry (Assumption A3), but row N (ROW as
origin, i.e. entering the US) and column N (ROW as destination, i.e. leaving
the US) are separate free parameters per nativity, not required to match --
`f_ROW_in` is the row-N cost (ROW→state), `f_ROW_out` is the column-N cost
(state→ROW). f_{N,N}=0 throughout (staying in ROW costs nothing).
"""
function Parameters(fixed::NamedTuple, fᵈ_ROW_in::Real, fᵈ_ROW_out::Real, fᶠ_ROW_in::Real, fᶠ_ROW_out::Real)

    (; β, r, θ, δ, ρ, μ_z, ξ_z, ξ_ω, νᵈ, νᶠ, A, fᵈ, fᶠ, L, N) = fixed

    fᵈ_full = zeros(N, N)
    fᶠ_full = zeros(N, N)
    fᵈ_full[1:L, 1:L] .= fᵈ
    fᶠ_full[1:L, 1:L] .= fᶠ
    fᵈ_full[1:L, N]   .= fᵈ_ROW_out
    fᵈ_full[N, 1:L]   .= fᵈ_ROW_in
    fᶠ_full[1:L, N]   .= fᶠ_ROW_out
    fᶠ_full[N, 1:L]   .= fᶠ_ROW_in

    return Parameters(β, r, θ, δ, ρ, μ_z, ξ_z, ξ_ω, νᵈ, νᶠ, A, fᵈ_full, fᶠ_full, N)

end

"""
Convenience wrapper: load every fixed input from disk, then build Parameters.
For one-off calls (SsSolve.jl's demos); a grid/threaded search should call
load_fixed_parameters() once and use the Parameters(fixed, ...) method above
instead, to avoid repeated disk I/O and JLD2's thread-unsafety.
"""
function Parameters(; fᵈ_ROW_in::Real, fᵈ_ROW_out::Real, fᶠ_ROW_in::Real, fᶠ_ROW_out::Real,
                     β::Real = 0.96)

    fixed = load_fixed_parameters(; β)
    return Parameters(fixed, fᵈ_ROW_in, fᵈ_ROW_out, fᶠ_ROW_in, fᶠ_ROW_out)

end

"""
Load 2015 US state labor supplies by nativity from StateAnalysisPreTfp.dta,
sorted by state FIPS code to match load_bilateral_costs()'s state ordering.
Covers the 51 US states only -- Rest of World's supply isn't in this file.
"""
function load_labor_supply()

    df = DataFrame(load(joinpath(data, "StateAnalysisPreTfp.dta")))
    df = df[df.year .== 2015, :]
    sort!(df, :statefip)

    return (Lᵈ = Float64.(df.Supply_Domestic), Lᶠ = Float64.(df.Supply_Foreign))

end

"""
Load 2015 US state wages by nativity from StateAnalysisPreTfp.dta, sorted by
state FIPS code to match load_bilateral_costs()'s state ordering.
"""
function load_wages()

    df = DataFrame(load(joinpath(data, "StateAnalysisPreTfp.dta")))
    df = df[df.year .== 2015, :]
    sort!(df, :statefip)

    return (Wᵈ = Float64.(df.Wage_Domestic), Wᶠ = Float64.(df.Wage_Foreign))

end

"""
Calibrate the output scale A (paper's Technology section) so that the
model's implied aggregate labor income, evaluated at the actual 2015 ACS
wages and labor supplies, equals the actual aggregate 2015 labor income.
Not identified by the migration data (which only pins down wage ratios), so
this is a units convention, not an estimating equation.
"""
function calibrate_scale(θ::Real, δ::Real, r::Real, ρ::Real, μ_z::Real, ξ_ω::Real, ξ_z::Real)

    (; Wᵈ, Wᶠ) = load_wages()
    (; Lᵈ, Lᶠ) = load_labor_supply()

    data_income = sum(Wᵈ .* Lᵈ .+ Wᶠ .* Lᶠ)

    model_income_at_A1 = 0.0
    for l in eachindex(Wᵈ)
        (; Z, λ, oneMinusλ) = TaskAggregates_LN(ρ, μ_z, ξ_ω, ξ_z, Wᵈ[l] / Wᶠ[l])
        L = LaborAggregate(λ, ρ, Lᶠ[l], Lᵈ[l], oneMinusλ)
        model_income_at_A1 += Z * L
    end
    model_income_at_A1 *= (1 - θ) * (θ / (r + δ))^(θ / (1 - θ))

    return data_income / model_income_at_A1

end

"""
Total 2015 world population by nativity, used to scale the stationary
labor-supply distribution: total_Lᵈ is everyone born in the US, wherever
they live; total_Lᶠ is everyone else in the world.

total_Lᵈ = 2015 US domestic-born labor supply (load_labor_supply(), summed
across the 51 states) + US-born people living abroad (MakePopulationConstants.do's
us_abroad_2015). total_Lᶠ = world population (that file's world_pop_2015)
minus total_Lᵈ.
"""
function load_total_population()

    Lᵈ_us = sum(load_labor_supply().Lᵈ)
    pop   = only(eachrow(DataFrame(load(joinpath(data, "PopulationConstants.dta")))))
    total_Lᵈ = Lᵈ_us + pop.us_abroad_2015
    total_Lᶠ = pop.world_pop_2015 - total_Lᵈ

    return (; total_Lᵈ, total_Lᶠ)

end

"""
Solve Vᵈ, Vᶠ, Lᵈ, Lᶠ, Πᵈ, Πᶠ given wages Wᵈ, Wᶠ (each length N) and the total
domestic/foreign population (US + Rest of World) used to scale the
stationary labor-supply distribution.
"""
function Solution(p::Parameters, Wᵈ::Vector{<:Real}, Wᶠ::Vector{<:Real};
                   total_Lᵈ::Real, total_Lᶠ::Real)
    Wᵈ, Wᶠ = Float64.(Wᵈ), Float64.(Wᶠ)
    Vᵈ, Πᵈ = solve_value_choiceprobs(Wᵈ, p.fᵈ, p.β, p.νᵈ)
    Vᶠ, Πᶠ = solve_value_choiceprobs(Wᶠ, p.fᶠ, p.β, p.νᶠ)
    Lᵈ = stationary_distribution(Πᵈ, total_Lᵈ)
    Lᶠ = stationary_distribution(Πᶠ, total_Lᶠ)
    return Solution(Vᵈ, Vᶠ, Wᵈ, Wᶠ, Lᵈ, Lᶠ, Πᵈ, Πᶠ)
end
