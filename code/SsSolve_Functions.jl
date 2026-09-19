using DataFrames, StatFiles, JLD2, LinearAlgebra

include("Estimation_Funcs.jl")

# %% Structs =================================================================

"""
Model parameters. fᵈ, fᶠ are N×N bilateral migration cost matrices, N = 51 US
states + Rest of World (the last row/column).
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
    fᵈ::Matrix{T}
    fᶠ::Matrix{T}
    N::Int
end

"""
Steady-state values, wages, and choice probabilities, by nativity. Wᵈ, Wᶠ are
inputs here, not yet solved for by market clearing.
"""
struct Solution{T <: Real}
    Vᵈ::Vector{T}
    Vᶠ::Vector{T}
    Wᵈ::Vector{T}
    Wᶠ::Vector{T}
    Πᵈ::Matrix{T}
    Πᶠ::Matrix{T}
end

# %% Solver Functions =========================================================

logsumexp(x) = (m = maximum(x); m + log(sum(exp.(x .- m))))
softmax(x)   = (e = exp.(x .- maximum(x)); e ./ sum(e))

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
Build Parameters. Every value except fᵈ_ROW, fᶠ_ROW (the domestic/foreign
migration cost between any US state and Rest of World, assumed common across
states) is loaded from an already-estimated source.
"""
function Parameters(; fᵈ_ROW::Real, fᶠ_ROW::Real, β::Real = 0.96)

    r  = 1 / β - 1
    θδ = load_theta_delta()
    θ, δ = θδ.theta_l[1], θδ.delta_l[1]
    ρ = load_aggsupply_estimate().ρ
    (; μ_z, ξ_ω, ξ_z) = load_cp_estimate()
    (; νᵈ, νᶠ) = load_scale_estimate()
    (; fᵈ, fᶠ) = load_bilateral_costs()

    L = size(fᵈ, 1)
    N = L + 1

    fᵈ_full = zeros(N, N)
    fᶠ_full = zeros(N, N)
    fᵈ_full[1:L, 1:L] .= fᵈ
    fᶠ_full[1:L, 1:L] .= fᶠ
    fᵈ_full[1:L, N]   .= fᵈ_ROW
    fᵈ_full[N, 1:L]   .= fᵈ_ROW
    fᶠ_full[1:L, N]   .= fᶠ_ROW
    fᶠ_full[N, 1:L]   .= fᶠ_ROW

    return Parameters(Float64(β), Float64(r), Float64(θ), Float64(δ), Float64(ρ),
                       Float64(μ_z), Float64(ξ_z), Float64(ξ_ω), Float64(νᵈ), Float64(νᶠ),
                       Matrix{Float64}(fᵈ_full), Matrix{Float64}(fᶠ_full), N)

end

"""
Solve Vᵈ, Vᶠ, Πᵈ, Πᶠ given wages Wᵈ, Wᶠ (each length N).
"""
function Solution(p::Parameters, Wᵈ::Vector{<:Real}, Wᶠ::Vector{<:Real})
    Wᵈ, Wᶠ = Float64.(Wᵈ), Float64.(Wᶠ)
    Vᵈ, Πᵈ = solve_value_choiceprobs(Wᵈ, p.fᵈ, p.β, p.νᵈ)
    Vᶠ, Πᶠ = solve_value_choiceprobs(Wᶠ, p.fᶠ, p.β, p.νᶠ)
    return Solution(Vᵈ, Vᶠ, Wᵈ, Wᶠ, Πᵈ, Πᶠ)
end
