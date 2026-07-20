"""
Use the expm1() function to write ω̄ᵖ - 1 where p is some power. This function is much
more stable computationally. Further, its derivatives can be smoothly approximated using the Taylor Series
for the exponential in neighborhoods where y → 0.
"""
function stable_expm1_ratio(y)
    ay = y isa ForwardDiff.Dual ? abs(ForwardDiff.value(y)) : abs(y)
    if ay < 1e-5
        return 1 + y/2 + y^2/6 + y^3/24 + y^4/120
    else
        return expm1(y) / y
    end
end

"""
Foreign-born integral (document's Integral 1): I_F = w^(-b)[D + m/(1-b)] where
D := m(1-m^(a-1))/(a-1), stabilized via stable_expm1_ratio at the a=1 knife edge.
"""
function I_F_μ(ρ, γ, μ, w)

    b = ρ / (1 - ρ)
    a = 1/γ + b
    m = μ * w
    lnm = log(m)
    y = (a - 1) * lnm

    D = -m * lnm * stable_expm1_ratio(y)

    return w^(-b) * (D + m / (1 - b))
end

"""
Domestic-born integral: I_D = [1 - m - D] / (1+γb), same D as I_F_μ.
"""
function I_D_μ(ρ, γ, μ, w)

    b = ρ / (1 - ρ)
    a = 1/γ + b
    m = μ * w
    lnm = log(m)
    y = (a - 1) * lnm

    D = -m * lnm * stable_expm1_ratio(y)

    return (1 - m - D) / (1 + γ * b)
end

"""
Task productivity aggregate Z and foreign-born task share λ, free-μ parameterization.
"""
function TaskAggregates_μ(ρ, γ, μ, w)

    b = ρ / (1 - ρ)
    F = I_F_μ(ρ, γ, μ, w)
    D = I_D_μ(ρ, γ, μ, w)

    Z = (F + D)^(1 / b)
    λ = F / (F + D)

    return (; Z, λ)
end

"""
CES aggregate of foreign and domestic labor supplies given the foreign-born task share λ.
"""
function LaborAggregate(λ, ρ, LF, LD)

    return (λ^(1 - ρ) * LF^ρ + (1 - λ)^(1 - ρ) * LD^ρ)^(1 / ρ)
end

"""
Sum of squared residuals for equation for production function (eqn 4.2 of the text)
    - Time FEs ι_t, first year normalized to 0, ι is the intercept
    - μ is the minimum of the Pareto
    - θ is the capital share
    - ρ controls the EOS between domestic and foreign

    x = [ι, θ, μ, γ, ρ, ιₜ_2,...,ιₜ_T];

    time_id maps each observation to an index into ιₜ (1..T).
    Returns Inf if 1 ≤ρ/(1-ρ) (violates the Pareto tail-moment condition, under which
    Z's defining integral fails to converge)
"""
function RSS(x::AbstractVector, time_id::AbstractVector{<:Integer}, T::Integer,
             Y, K, wD, wF, LF, LD)

    ι, θ, μ, γ, ρ = x[1], x[2], x[3], x[4], x[5]
    ιₜ = vcat(zero(ι), x[6:4+T])
    b = ρ / (1 - ρ)
    b < 1 || return oftype(ρ, Inf)

    w = wD ./ wF
    ta = TaskAggregates_μ.(ρ, γ, μ, w)
    Z = getproperty.(ta, :Z)
    λ = getproperty.(ta, :λ)
    L = LaborAggregate.(λ, ρ, LF, LD)

    ŷ = ι .+ ιₜ[time_id] .+ θ .* log.(K) .+ (1 - θ) .* (log.(Z) .+ log.(L))
    res = log.(Y) .- ŷ

    return sum(abs2, res)
end

"""
Estimate the production function via Fminbox(LBFGS) directly on x = [ι, θ, μ, γ, ρ,
ιₜ_2,...,ιₜ_T], with time fixed effects.
"""
function EstimateProdFunc(x0::AbstractVector, lb::AbstractVector, ub::AbstractVector,
                           time_id::AbstractVector{<:Integer}, T::Integer,
                           Y, K, wD, wF, LF, LD; iterations = 1000, show_trace = false, g_tol = 1e-8)

    obj(x) = RSS(x, time_id, T, Y, K, wD, wF, LF, LD)
    opts = Optim.Options(iterations = iterations, show_trace = show_trace, g_tol = g_tol)

    result = Optim.optimize(obj, lb, ub, x0, Optim.Fminbox(Optim.LBFGS()), opts; autodiff = :forward)

    xhat = Optim.minimizer(result)
    ι, θ, μ, γ, ρ = xhat[1], xhat[2], xhat[3], xhat[4], xhat[5]
    ιₜ = vcat(zero(ι), xhat[6:4+T])

    return (; result, ι, θ, μ, γ, ρ, ιₜ)
end