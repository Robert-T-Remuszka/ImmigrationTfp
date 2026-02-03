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
    ψ::T1                                       # Parameterization for χ: χ > γᵈ/(γᵈ - γᶠ) ⟺ χ > (γᶠ + Δ)/Δ ⟺ χ = γᶠ/Δ + 1 + exp(ψ)
    ι::T1                                       # Intercept of production function
    ιₛ::Vector{T1}                               # State fixed effects
    ιₜ::Vector{T1}                               # Time fixed effects
    
end

"""
Constructor for the AuxParameters type.
"""
function AuxParameters(;
    ρ::T1  = 0.50,
    θ::T1  = 0.50,
    γᶠ::T1 = 1.,
    Δ::T1  = 1.,
    Γ::T1  = -1.,
    αᶠ::T1 = 2.,
    αᵈ::T1 = 4.,
    ψ::T1  = 2.,
    ι::T1  = 5.,
    df::DataFrame   = StateAnalysis,
    N::Int          = length(unique(df[:, :statefip])),
    T::Int          = length(unique(df[:, :year])),
    ιₛ::Vector{T1} = zeros(N),
    ιₜ::Vector{T1} = zeros(T)
    ) where{T1 <: Real}                
    
    return AuxParameters{T1}(ρ, θ, γᶠ, Δ, Γ, αᶠ, αᵈ, ψ, ι, ιₛ, ιₜ)

end
 
"""
Compute parts of the reduced form production function for given parameters p.
"""
function ComputeReduced(wᶠ::T1, wᵈ::T1, F::T1, D::T1; p::AuxParameters) where {T1 <: Real}
    
    (; ρ, γᶠ, Δ, Γ, αᶠ, αᵈ, ψ) = p

    b, γᵈ  = ρ / (1 - ρ), γᶠ + Δ
    bᶠ, bᵈ, b_gamma = γᶠ * b, γᵈ * b, Γ * b
    χ      = γᶠ / Δ + 1 + exp(ψ)
    w = wᵈ / wᶠ
    α = αᶠ / αᵈ
    Z = max(
        (1/bᶠ) * ((α * w / exp(Γ))^(bᶠ/Δ) * (Δ * χ)/(Δ * χ - γᵈ) - χ/(χ - 1)) + 
        exp(b_gamma) * (1/bᵈ) * (exp(bᵈ) - (α * w / exp(Γ))^(bᵈ / Δ) * (Δ * χ)/(Δ * χ - γᵈ))
        , 0.01
        )^(1 / b)
    λ = clamp( (1/bᶠ) * ((α * w / exp(Γ))^(bᶠ/Δ) * (Δ * χ)/(Δ * χ - γᵈ) - χ/(χ - 1)) / Z^b, 0. , 1.)
    L = (λ^(1 - ρ) * (αᶠ * F)^ρ + (1 - λ)^(1 - ρ) * (αᵈ * D)^ρ)^(1/ρ)

    return (; Z, λ, L)

end

"""
Compute residual sum of suares for a set of parameters p.
"""
function RSS(x::Vector{T1}; df::DataFrame = StateAnalysis) where {T1 <: Real}
    
    # Ensure the data is sorted correctly - Matters for fixed effects
    df_sort = @chain df begin
        @arrange(statefip, year)
    end

    N, T = length(unique(df_sort[:, :statefip])), length(unique(df_sort[:, :year]))

    # Package the input vector into instance of AuxParameters. A leading zero is placed in front fixed effect vecors
    # since I include an intercept.
    p = AuxParameters(
        ρ   = x[1],
        θ   = x[2], 
        γᶠ  = x[3], 
        Δ   = x[4], 
        Γ   = x[5], 
        αᶠ  = x[6], 
        αᵈ  = x[7], 
        ψ   = x[8],
        ι   = x[9], 
        ιₛ = vcat(0., x[10 : 8 + N]),
        ιₜ = vcat(0., x[9 + N : end]),
        df  = df_sort
    )

        
    # Unpack some necessary parameters
    (; ι, θ, ιₛ, ιₜ) = p
    state_indices = repeat(1:N, inner=T)  # [1,1,...,1, 2,2,...,2, ..., N,N,...,N]
    time_indices = repeat(1:T, outer=N)   # [1,2,...,T, 1,2,...,T, ..., 1,2,...,T]
    state_fe = ιₛ[state_indices]
    time_fe = ιₜ[time_indices]

    # Fetch some data for estimation
    Y, K         = df_sort[:, :GDP], df_sort[:, :CapStock]
    wᶠ, wᵈ, F, D = df_sort[:, :Wage_Foreign], df_sort[:, :Wage_Domestic], df_sort[:, :Supply_Foreign], df_sort[:, :Supply_Domestic]

    # Broadcastisting ComputeReduced returns a vector of named tuples, so we need to broadcast the getproperty function
    # getproperty is the basic function for which '.' notation is a shorthand, but . is also a shorthand for broadcasting!
    Z, L = getproperty.(ComputeReduced.(wᶠ, wᵈ, F, D; p = p), :Z), getproperty.(ComputeReduced.(wᶠ, wᵈ, F, D; p = p), :L)
    res  = log.(Y) - (ι .+  state_fe .+  time_fe .+ θ * log.(K) + (1 - θ) * log.(Z .* L))

    return dot(res, res)

end

"""
Estimate production function by minimizing residual sum of squares. Initial guess is x0.
"""
function EstimateProdFunc(x0::Vector{T1}; df::DataFrame = StateAnalysis) where {T1 <: Real}

    obj(x) = RSS(x; df = df)
    opts   = Optim.Options(
        show_trace = true, 
        show_every = 50,
        f_tol = 1e-8,        # Function tolerance
        x_tol = 1e-8,        # Parameter tolerance
        g_tol = 1e-6,  
        iterations = 100000   # Allow many more iterations
    )
    N, T   = length(unique(df[:, :statefip])), length(unique(df[:, :year]))
    
    lb, ub = zeros(length(x0)), zeros(length(x0))
    lb[1] = 0.01          # ρ > 0
    ub[1] = 0.99          # ρ < 1
    lb[2] = 1e-2          # θ > 0
    ub[2] = 1. - 1e-2     # θ < 1
    lb[3] = 1e-1          # γᶠ > 0
    ub[3] = Inf
    lb[4] = 1e-1          # Δ > 0
    ub[4] = Inf
    lb[5] = -Inf          # Γ unbounded
    ub[5] =  Inf
    lb[6] =  0.05         # αᶠ > 0
    ub[6] =  Inf
    lb[7] =  0.05         # αᵈ > 0
    ub[7] =  Inf
    lb[8] = -Inf          # ψ unbounded
    ub[8] =  Inf
    lb[9] = -Inf          # ι unbounded
    ub[9] =  Inf
    lb[10 : 8 + N] .= -Inf  # ιₛ unbounded
    ub[10 : 8 + N] .=  Inf
    lb[9 + N : end] .= -Inf # ιₜ unbounded
    ub[9 + N : end] .=  Inf

    result = optimize(obj, lb, ub, x0, Fminbox(NelderMead()), opts)
    MSE    =  result.minimum / (N * T)

    return (; result.minimizer, MSE)

end