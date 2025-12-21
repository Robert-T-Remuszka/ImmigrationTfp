"""
Parameters of the auxilliary model for TFP. I call this an auxilliary model because
it could be repeated on model-simulated data.
"""
struct AuxParameters{T1 <: Real}

    ρ::T1                                       # CES parameter
    θ::T1                                       # Capital share
    γᶠ::T1                                      # Comp advantage - foreign
    Δ::T1                                       # γᵈ = γᶠ + Δ
    Γ::T1                                       # Comp advantage - level shifter
    αᶠ::T1                                      # Absolute advantage - foreign
    αᵈ::T1                                      # Absolute advantage - domestic
    ι::T1                                       # Intercept of production function
    
end

"""
Constructor for the AuxParameters type.
"""
function AuxParameters(;
    ρ::T1  = 0.25,
    θ::T1  = 0.50,
    γᶠ::T1 = 1.,
    Δ::T1  = 2.,
    Γ::T1  = -1.,
    αᶠ::T1 = 3.,
    αᵈ::T1 = 3.,
    ι::T1  = 5.) where{T1 <: Real}                
    
    return AuxParameters{T1}(ρ, θ, γᶠ, Δ, Γ, αᶠ, αᵈ, ι)

end
 
"""
Compute parts of the reduced form production function for given parameters p.
"""
function ComputeReduced(wᶠ::T1, wᵈ::T1, F::T1, D::T1; p::AuxParameters) where {T1 <: Real}
    
    (; ρ, γᶠ, Δ, Γ, αᶠ, αᵈ) = p

    b, γᵈ  = ρ / (1 - ρ), γᶠ + Δ
    bᶠ, bᵈ, b_gamma = γᶠ * b, γᵈ * b, Γ * b
    w = wᵈ ./ wᶠ
    α = αᶠ  / αᵈ
    𝒯 = clamp((log(w * α) - Γ) / Δ, 0., 1.)
    Z = max((1 / bᶠ) * (exp(bᶠ * 𝒯) - 1) + (1 / bᵈ) * exp(b_gamma) * (exp(bᵈ) - exp(bᵈ * 𝒯)), 0.)^(1 / b)
    λ = clamp((1 / bᶠ) * (exp(bᶠ * 𝒯) - 1) / Z^b, 0. , 1.)
    L = (λ^(1 - ρ) * (αᶠ * F)^ρ + (1 - λ)^(1 - ρ) * (αᵈ * D)^ρ)^(1/ρ)

    return (; 𝒯, Z, λ, L)

end

"""
Compute residual sum of suares for a set of parameters p.
"""
function RSS(x::Vector{T1}; df::DataFrame = StateAnalysis) where {T1 <: Real}
    
    p = AuxParameters(ρ = x[1], θ = x[2], γᶠ = x[3], Δ = x[4], Γ = x[5], αᶠ = x[6], αᵈ = x[7], ι = x[8])
    (; ι, θ) = p

    Y, K         = df[:, :GDP] , df[:, :CapStock]
    wᶠ, wᵈ, F, D = df[:, :Wage_Foreign], df[:, :Wage_Domestic], df[:, :Supply_Foreign], df[:, :Supply_Domestic]

    # Broadcastisting ComputeReduced returns a vector of named tuples, so we need to broadcast the getproperty function
    # getproperty is the basic function for which '.' notation is a shorthand, but . is also a shorthand for broadcasting!
    Z, L = getproperty.(ComputeReduced.(wᶠ, wᵈ, F, D; p = p), :Z), getproperty.(ComputeReduced.(wᶠ, wᵈ, F, D; p = p), :L)
    res  = log.(Y) - (ι .+ θ * log.(K) + (1 - θ) * log.(Z .* L))

    return dot(res, res)

end

"""
Estiamte production function by minimizing residual sum of squares. Initial guess is x0.
"""
function EstimateProdFunc(x0::Vector{T1}; df::DataFrame = StateAnalysis) where {T1 <: Real}

    obj(x) = RSS(x; df = df)
    opts   = Optim.Options(show_trace = true, g_tol = 1e-6, show_every = 10)
    N, T   = length(unique(df[:, :statefip])), length(unique(df[:, :year]))
    
    lb, ub = zeros(length(x0)), zeros(length(x0))
    lb[1] = 1e-2          # ρ > 0
    ub[1] = 1. - 1e-2     # ρ < 1
    lb[2] = 1e-2          # θ > 0
    ub[2] = 1. - 1e-2     # θ < 1
    lb[3] = 1e-4          # γᶠ > 0
    ub[3] = Inf
    lb[4] = 1e-4          # Δ > 0
    ub[4] = Inf
    lb[5] = -Inf          # Γ unbounded
    ub[5] =  Inf
    lb[6] =  1e-4         # αᶠ > 0
    ub[6] =  Inf
    lb[7] =  1e-4         # αᵈ > 0
    ub[7] =  Inf
    lb[8] = -Inf          # ι unbounded
    ub[8] =  Inf

    result = optimize(obj, lb, ub, x0, Fminbox(NelderMead()), opts)
    MSE    =  result.minimum / (N * T)

    return (; result.minimizer, MSE)

end