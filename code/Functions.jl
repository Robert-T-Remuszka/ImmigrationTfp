"""
Parameters of the auxilliary model for TFP. I call this an auxilliary model because
it could be repeated on model-simulated data.
"""
struct AuxParameters{T1 <: Real}

    ρ::T1                                       # CES parameter
    θ::T1                                       # Capital share
    ζᶠ::T1                                      # Comp advantage - foreign
    ζᵈ::T1                                      # Comp advantage - domestic
    αᶠ::T1                                      # Absolute advantage - foreign, normalize αᵈ = 1.
    ξ::Vector{T1}                               # Time fixed effects - national task measure control
    χ::Vector{T1}                               # State fixed effects
    Inter::T1                                   # Intercept
    
end

"""
Constructor for the AuxParameters type.
    - df should be sorted, state then year.
    - It is assumed that domestic have comp advantage in certain tasks. This assumption is
    reflected in the parameterization ζᵈ = ζᶠ + Δ
    - I will also normalize αᶠ < αᵈ = 1. The idea is that it is only the ratio of the alphas
    that matters (task cutoff and CRS L). This will be paramterized by αᶠ = 1 / (1 + exp(δ)) where
    δ is a real number.
"""
function AuxParametersConstructor(;
    ρ::T1 = 0.5,
    θ::T1 = 0.3,
    ζᶠ::T1 = 1.5,
    Δ::T1 = 5.,
    δ::T1 = 1., 
    df::DataFrame = StateAnalysis,
    T::Int = length(unique(df[:,:year])),
    N::Int = length(unique(df[:,:statefip])),
    ξ::Vector{T1} = zeros(T - 1),
    χ::Vector{T1} = zeros(N -1),
    Inter::T1 = 0.
    ) where{T1 <: Real}                
    
    return AuxParameters{T1}(ρ, θ, ζᶠ, ζᶠ + Δ, 1 / (1 + exp(δ)), ξ, χ, Inter)

end

"""
Calculate the cutoff task. The notation, s̄ comes from a change of variables.
"""
function s̄(p::AuxParameters; df::DataFrame = StateAnalysis)

    (; ζᶠ, ζᵈ, αᶠ, ρ) = p
    
    # Unpack wages
    b = ρ / (1 - ρ)
    #ζ = ((1 + ζᶠ * b)  / (1 + ζᵈ * b))^(1/b)
    wᵈ = df[:, :Wage_Domestic] ./  mean(df[:, :Wage_Domestic])
    wᶠ = df[:, :Wage_Foreign] ./ mean(df[:, :Wage_Foreign])
    w = wᵈ ./ wᶠ
    w_level = w .* (mean(df[:, :Wage_Domestic]) ./ mean(df[:, :Wage_Foreign]))
    s̄ = (w * αᶠ).^(1/(ζᵈ - ζᶠ))
    s̄_level = (w_level * αᶠ).^(1/(ζᵈ - ζᶠ))
    
    return s̄, s̄_level

end

"""
The residual function. 
"""
function Residual(p::AuxParameters; df::DataFrame = StateAnalysis)
    
    # Unpacking
    (; ρ, θ, ζᶠ, ζᵈ, αᶠ, ξ, χ, Inter) = p
    K = df[:, :CapStock] ./ mean(df[:, :CapStock])
    F = df[:, :Supply_Foreign] ./ mean(df[:, :Supply_Foreign])
    D = df[:, :Supply_Domestic] ./ mean(df[:, :Supply_Domestic])
    D_level = df[:, :Supply_Domestic]
    F_level = df[:, :Supply_Foreign]
    Y = df[:, :GDP] ./ mean(df[:, :GDP])
    T = length(unique(df[:, :year]))
    N = length(unique(df[:,:statefip]))
    TFE_Mat = repeat(Matrix{Float64}(I, T, T), N)
    ξ_conform = vcat(0., ξ)
    χ_conform = vcat(0., χ)

    # Set up the state fixed effects matrix
    SFE_Mat = Matrix{Float64}(undef, 0, N)
    for c in 1:N
        SFE_Mat = vcat(SFE_Mat, [j == c ? 1. : 0. for i in 1:T, j in 1:N])
    end

    # Calculate parts of the production function
    b = ρ/(1 - ρ)
    s = s̄(p; df = df)
    Z = max.(s[1].^(1 + ζᶠ * b) ./ (1 + ζᶠ * b)  .+ (1 .- s[1].^(1 + ζᵈ * b)) ./ (1 + ζᵈ * b), 0.)
    Z_level = max.(s[2].^(1 + ζᶠ * b) ./ (1 + ζᶠ * b)  .+ (1 .- s[2].^(1 + ζᵈ * b)) ./ (1 + ζᵈ * b), 0.)
    λ = clamp.(s[1].^(1 + ζᶠ * b) ./ (1 + ζᶠ * b) ./ Z, 0. , 1.)
    λ_level = clamp.(s[2].^(1 + ζᶠ * b) ./ (1 + ζᶠ * b) ./ Z_level, 0. , 1.)
    L = (λ.^(1 - ρ) .* (αᶠ * F).^ρ + (1 .- λ).^(1 - ρ) .* D.^ρ).^(1/ρ)
    L_level = (λ_level.^(1 - ρ) .* (αᶠ * F_level).^ρ + (1 .- λ_level).^(1 - ρ) .* D_level.^ρ).^(1/ρ)
    

    return log.(Y) - (Inter .+ TFE_Mat * ξ_conform  + SFE_Mat * χ_conform + θ * log.(K) + (1 - θ)/b * log.(Z) + (1 - θ) * log.(L)), Z_level, L_level

end

"""
The sum of squared errors for production function.
"""
function SSE(x::Vector{T1}; df::DataFrame = StateAnalysis) where{T1 <: Real}

    ρ, θ, ζᶠ, Δ, δ, ξ, χ, Inter = x[1], x[2], x[3], x[4], x[5], x[6 : 4 + T], x[5 + T : 3 + T + N], x[end]

    p = AuxParametersConstructor(ρ = ρ, θ = θ, ζᶠ = ζᶠ, Δ = Δ, δ = δ, ξ = ξ, χ = χ, Inter = Inter)
    vals = Residual(p; df = df)

    return dot(vals[1], vals[1])

end

"""
Set up and call the box constrained optimization
"""
function EstimateProduction(x0::Vector{T1}; df::DataFrame = StateAnalysis) where{T1 <: Real}

    N, T = length(unique(df[:,:statefip])) ,length(unique(df[:,:year]))
    options = Optim.Options(outer_iterations = 10000, iterations = 10000, show_trace = true, show_every = 10, g_tol = 1e-6);

    # Set up the constraints
    lb = zeros(5 + T - 1 + N - 1 + 1);                         # Initialize all bounds first
    ub = Inf * ones(5 + T - 1 + N - 1 + 1);                    # Initialize all bounds first

    # Implment parameter-specific bounds
    ub[1:2] .= 1.;                                           # ρ, θ < 1
    ub[3] = 12.
    ub[4] = 20.
    lb[5] = -2.;                                             # Don't want δ too low or s̄ > 1
    ub[5] = 20.
    lb[6: 4 + T] .= -Inf;                                    # Unrestricted Time FEs
    lb[5 + T : 3 + T + N] .= -Inf;                           # Unrestricted State FEs
    lb[end] = -Inf;                                          # The intercept is unrestricted
    
    # Objective and gradient
    f(x) = SSE(x; df = df)
    function g!(G, x)
        ForwardDiff.gradient!(G, x -> SSE(x; df=df), x)
    end

    # Call the optimizer
    res = optimize(f, g!, lb, ub, x0, Fminbox(LBFGS()), options)
    MSE = res.minimum / (N * T)

    return res.minimizer, MSE

end