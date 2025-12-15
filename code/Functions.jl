"""
Parameters of the auxilliary model for TFP. I call this an auxilliary model because
it could be repeated on model-simulated data.
"""
struct AuxParameters{T1 <: Real}

    ρ::T1                                       # CES parameter
    θ::T1                                       # Capital share
    ζᶠ::T1                                      # Comp advantage - foreign
    ζᵈ::T1                                      # Comp advantage - domestic
    αᶠ::T1                                      # Absolute advantage - foreign
    αᵈ::T1                                      # Absolute advantage - domestic
    μ::T1                                       # Task measure

end

"""
Constructor for the AuxParameters type.
    - df should be sorted, state then year.
    - It is assumed that domestic have comp advantage in certain tasks. This assumption is
    reflected in the parameterization ζᵈ = ζᶠ + Δ
"""
function AuxParametersConstructor(;
    ρ::T1 = 0.5,
    θ::T1 = 0.3,
    ζᶠ::T1 = 1.5,
    Δ::T1 = 5.,
    αᶠ::T1 = 1.,
    αᵈ::T1 = 2.,
    μ::T1 = 10.
    ) where{T1 <: Real}                
    
    return AuxParameters{T1}(ρ, θ, ζᶠ, ζᶠ + Δ, αᶠ, αᵈ, μ)

end

"""
Calculate the cutoff task.
"""
function T(p::AuxParameters; df::DataFrame = StateAnalysis)

    (; ζᶠ, ζᵈ, αᶠ, αᵈ, μ) = p
    
    wᵈ, wᵈ_lvl = df[:, :Wage_Domestic] ./  std(df[:, :Wage_Domestic]), df[:, :Wage_Domestic]
    wᶠ, wᶠ_lvl = df[:, :Wage_Foreign]  ./  std(df[:, :Wage_Foreign]),  df[:, :Wage_Foreign]
    w, w_lvl   = wᵈ ./ wᶠ, wᵈ_lvl ./ wᶠ_lvl
    α          = αᶠ ./ αᵈ
    τ, τ_lvl   = (w * α).^(1/(ζᵈ - ζᶠ)), (w_lvl * α).^(1/(ζᵈ - ζᶠ))
    
    return clamp.(τ, 0. , μ), clamp.(τ_lvl, 0., μ)

end

"""
The residual function. 
"""
function Residual(p::AuxParameters; df::DataFrame = StateAnalysis)
    
    # Unpacking
    (; ρ, θ, ζᶠ, ζᵈ, αᶠ, αᵈ, μ) = p
    K, K_lvl = df[:, :CapStock] ./ std(df[:, :CapStock]),               df[:, :CapStock]
    F, F_lvl = df[:, :Supply_Foreign] ./ std(df[:, :Supply_Foreign]),   df[:, :Supply_Foreign]
    D, D_lvl = df[:, :Supply_Domestic] ./ std(df[:, :Supply_Domestic]), df[:, :Supply_Domestic]
    Y = df[:, :GDP] ./ std(df[:, :GDP])

    # Calculate parts of the production function
    b        = ρ/(1 - ρ)
    τ, τ_lvl = T(p; df = df)
    Z, Z_lvl = (τ.^(1 + ζᶠ * b) ./ (1 + ζᶠ * b) + (μ^(1 + ζᵈ * b) .- τ.^(1 + ζᵈ * b)) ./ (1 + ζᵈ * b)).^(1/b), (τ_lvl.^(1 + ζᶠ * b) ./ (1 + ζᶠ * b) + (μ^(1 + ζᵈ * b) .- τ_lvl.^(1 + ζᵈ * b)) ./ (1 + ζᵈ * b)).^(1/b)
    λ, λ_lvl = clamp.((τ.^(1 + ζᶠ * b) ./ (1 + ζᶠ * b)) ./ Z.^b, 0. , 1.), clamp.((τ_lvl.^(1 + ζᶠ * b) ./ (1 + ζᶠ * b)) ./ Z_lvl.^b, 0. , 1.)
    L, L_lvl = (λ.^(1 - ρ) .* (αᶠ * F).^ρ + (1 .- λ).^(1 - ρ) .* D.^ρ).^(1/ρ), (λ_lvl.^(1 - ρ) .* (αᶠ * F_lvl).^ρ + (1 .- λ_lvl).^(1 - ρ) .* D_lvl.^ρ).^(1/ρ)
    

    return log.(Y) - (θ * log.(K) + (1 - θ) * log.(Z .* L)), Z_lvl, L_lvl

end

"""
The sum of squared errors for production function.
"""
function SSE(x::Vector{T1}; df::DataFrame = StateAnalysis) where{T1 <: Real}

    ρ, θ, ζᶠ, Δ, αᶠ, αᵈ, μ = x[1], x[2], x[3], x[4], x[5], x[6], x[7]

    p = AuxParametersConstructor(ρ = ρ, θ = θ, ζᶠ = ζᶠ, Δ = Δ, αᶠ = αᶠ, αᵈ = αᵈ, μ = μ)
    vals = Residual(p; df = df)

    return dot(vals[1], vals[1])

end

"""
Set up and call the box constrained optimization
"""
function EstimateProduction(x0::Vector{T1}; df::DataFrame = StateAnalysis) where{T1 <: Real}

    # Initialize
    NS, NT = length(unique(df[:,:statefip])) ,length(unique(df[:,:year]))
    options = Optim.Options(outer_iterations = 10000, iterations = 10000, show_trace = true, show_every = 10, g_tol = 1e-6);
    lb = zeros(7);                         
    ub = zeros(7);                        

    # Implment parameter-specific bounds
    lb[1:2] .= 0.                            # ρ, θ > 0
    ub[1:2] .= 1.                            # ρ, θ < 1
    lb[3]    = 0.                            # ζᶠ > 0
    ub[3]    = 15.                                              
    lb[4]    = 0.                            # Δ > 0
    ub[4]    = 15.
    lb[5]    = 0.                           # Loosing nonnegativity on δ (if its -Inf the s̄ > 1 though, so not too loose)
    ub[5]    = 15.
    lb[6]    = 0.
    ub[6]    = 15.
    lb[7]    = 0.                                 
    ub[7]    = 15.
    
    # Objective and gradient
    f(x) = SSE(x; df = df)
    function g!(G, x)
        ForwardDiff.gradient!(G, x -> SSE(x; df=df), x)
    end

    # Call the optimizer
    res = optimize(f, lb, ub, x0, Fminbox(LBFGS()), options)
    MSE = res.minimum / (NS * NT)

    return res.minimizer, MSE

end