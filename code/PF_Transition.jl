#================================================================
                        TYPE DECLARATIONS
================================================================#
struct Parameters{T1 <: Real, T2 <: Integer}

    β::T1                                       # HH discount rate
    r::T1                                       # Capital rental rate
    δ::T1                                       # Capital depreciation rate
    ρ::T1                                       # CES parameter
    θ::T1                                       # Capital share
    γᶠ::T1                                      # Comp advantage - foreign
    γᵈ::T1                                      # Recall that γᵈ = γᶠ + Δ from the prelimenary estimation
    Γ::T1                                       # Comp advantage - level shifter
    αᶠ::T1                                      # Absolute advantage - foreign
    αᵈ::T1                                      # Absolute advantage - domestic
    ψ::T1                                       # Autoregressive coefficient
    M::Matrix{T1}                               # Quota policy - rows are origin, cols are destination
    νᵈ::T1                                      # Gumbel shape - domestic
    νᶠ::T1                                      # Gumbel shape - foreign

    # Data
    wᵈ_row::T1                                  # Domestic wage earned in 'rest of world' - thousands
    wᶠ_row::T1                                  # Foreign wage earned in 'rest of world' - thousands
    N::T2                                       # Number of regions; convention is that the N-th region is 'rest of the world'
    Πᵈ_pre::Matrix{T1}                          # Pre-period choice probabilities, domestic (rows origin, cols desintation)
    Πᶠ_pre::Matrix{T1}                          # Pre-period choice probabilities, foreign
    Lᵈ_tot::T1                                  # Total domestic population - billions
    Lᶠ_tot::T1                                  # Total foreign population  - billions

    T::T2                                       # Length of transition

end

struct TransitSoln{T1 <: Real}

    Wᵈ::Matrix{T1}                              # Domestic wages, rows are location and t is time
    Wᶠ::Matrix{T1}                              # Foreign wages
    Lᵈ::Matrix{T1}                              # Domestic labor supplies
    Lᶠ::Matrix{T1}                              # Foreign labor supplies
    U̇ᵈ::Matrix{T1}                              # Domestic value differences
    U̇ᶠ::Matrix{T1}                              # Foreign value changes
    Πᵈ::Array{T1}                               # Bilateral choice probabilities, domestic (rows origin, cols destination, height period)
    Πᶠ::Array{T1}                               # Bilateral choice probabilities, foreign

end

#================================================================
                        CONSTRUCTOR FUNCTIONS
================================================================#
function Parameters(p::AuxParameters = p_star;
    N::T1 = 3,
    β::T2  = 0.99,
    r::T2  = 1/β - 1,
    δ::T2  = 0.10,
    ρ::T2  = p.ρ,
    θ::T2  = p.θ,
    γᶠ::T2 = p.γᶠ,
    γᵈ::T2 = p.γᶠ + p.Δ,
    Γ::T2  = p.Γ,
    αᶠ::T2 = p.αᶠ,
    αᵈ::T2 = p.αᵈ,
    ψ::T2  = 0.45,
    National_Quota::T2 = (66000 + 85000) / 1e+9, # H-2B Statutory quota + H-1B Statutory with advanced degree exemption
    wᵈ_row::T2         = ( # See IRS Statistics of Income, deflated to 2009 real values using GDP deflator
      110297 / 1.24 # 2021
    + 119341 / 1.11 # 2016
    + 111807 / 1.03 # 2011
    + 111343 / 0.95 # 2006
    + 81379  / 0.84 # 2001
    + 60077  / 0.77 # 1996
    + 44090  / 0.70 # 1991
    ) / 7e+3,
    wᶠ_row::T2         = ( # See World Bank GNI and Population Statistics
      9328 / 1.24  # 2021
    + 8131 / 1.11  # 2016
    + 7842 / 1.03  # 2011
    + 5781 / 0.95  # 2006
    + 3869 / 0.84  # 2001
    + 4085 / 0.77  # 1996
    + 3391 / 0.70  # 1991
    ) / 7e+3,
    Lᵈ_tot::T2 = 0.3,
    Lᶠ_tot::T2 = 8.,
    νᵈ::T2 = 4.5,
    νᶠ::T2 = 4.5,
    T::T1 = 10
    ) where{T1 <: Integer, T2 <: Real}

    # make up some data for the pre-period migration shares - The rows are origin and cols are destination; this will come from data later on
    Πᵈ_pre = fill(1/N, N, N)
    Πᶠ_pre = fill(1/N, N, N)

    # Status-quo statutory quota policy - rows are origin, cols are destination
    M = fill(Inf64, N, N)
    M[N, 1:N-1] .= National_Quota / (N - 1)

    return Parameters{T2, T1}(β, r, δ, ρ, θ, γᶠ, γᵈ, Γ, αᶠ, αᵈ, ψ, M, νᵈ, νᶠ, wᵈ_row, wᶠ_row, N, Πᵈ_pre, Πᶠ_pre, Lᵈ_tot, Lᶠ_tot, T)

end

function TransitSoln(p::Parameters = Params1, Pct_domestic_abroad::T1 = 0.02, Pct_foreign_in_US::T1 = 0.007) where{T1 <: Real}

    (; N, Lᵈ_tot, Lᶠ_tot, wᵈ_row, wᶠ_row, Πᵈ_pre, Πᶠ_pre, T) = p
    
    # Initialize wages
    Wᵈ, Wᶠ = ones(N,T), ones(N, T)
    Wᵈ[N, :] .= wᵈ_row
    Wᶠ[N, :] .= wᶠ_row

    # Initialize choice probabilities
    Πᵈ, Πᶠ = fill(1/N, (N, N, T - 1)), fill(1/N, (N, N, T - 1))

    # Initialize supplies
    Lᵈ_pre = vcat(fill((1 - Pct_domestic_abroad) * Lᵈ_tot / (N - 1), N - 1), Pct_domestic_abroad * Lᵈ_tot) # These will be taken from the data, the pct parameters already are
    Lᶠ_pre = vcat(fill(Pct_foreign_in_US * Lᶠ_tot / (N - 1), N - 1), (1 - Pct_foreign_in_US) * Lᶠ_tot)
    Lᵈ, Lᶠ = zeros(N, T), zeros(N, T)
    for t in 1:T

        for l in 1:N
            
            for lp in 1:N

                Lᵈ[l,t] += t == 1 ? Πᵈ_pre[lp,l] * Lᵈ_pre[lp] : Πᵈ[lp,l,t - 1] * Lᵈ[lp, t - 1]
                Lᶠ[l,t] += t == 1 ? Πᶠ_pre[lp,l] * Lᶠ_pre[lp] : Πᶠ[lp,l,t - 1] * Lᶠ[lp, t - 1]

            end
            
        end

    end

    # Initialize value changes - the endoenous parts are the first T-1 columns, column T always stays at 1 to enforce the boundary condition
    U̇ᵈ, U̇ᶠ = ones(N, T), ones(N, T)

    return TransitSoln{T1}(Wᵈ, Wᶠ, Lᵈ, Lᶠ, U̇ᵈ, U̇ᶠ, Πᵈ, Πᶠ)

end

#================================================================
                        SOLVER FUNCTIONS
================================================================#
"""
Relative wage condition & resource feasability for a location at each time. Returns a 2 × T matrix of residuals.

Inputs:
    x := Relative wage
"""
function StaticEq(Lᵈ::Vector{T1}, Lᶠ::Vector{T1}, x::Vector{T1}, Wᶠ::Vector{T1}; p::Parameters = Params1) where{T1 <: Real}

    (; ρ, θ, r, δ, γᶠ, γᵈ, αᶠ, αᵈ, Γ) = p
    Wᵈ = x .* Wᶠ

    # Precompute the integral expressions
    Exp1         = ((γᵈ - γᶠ) *  (1 - ρ)) / ((γᵈ - γᶠ) * (1 - ρ) - ρ * γᵈ)                         # Eω^(ρ/(1-ρ) * γᵈ / (γᵈ - γᶠ))
    Exp2         = (1 - ρ) / (1 - 2 * ρ)                                                           # Eω^(ρ / (1-ρ))
    No_Arb_Stuff = (((αᶠ .* Wᵈ) ./ (αᵈ .* Wᶠ)) .* exp(-Γ)).^((ρ * γᵈ) / ((1 - ρ) * (γᵈ - γᶠ)))         # The coefficient coming from 𝒯

    # Compute parts of the production function
    Foreign_Int  = max.((1 - ρ) / (γᶠ * ρ) .* (No_Arb_Stuff .* Exp1 .- Exp2), 1e-4)
    Domestic_Int = max.(exp((ρ / (1 - ρ)) * Γ) * ((1 - ρ) / (ρ * γᵈ)) .* (exp((ρ / (1-ρ)) * γᵈ) .-  No_Arb_Stuff .* Exp1), 1e-4)
    Z, λ         = max.((Domestic_Int + Foreign_Int).^((1-ρ)/ρ), 1e-4), clamp.(Foreign_Int ./ (Foreign_Int + Domestic_Int), 1e-4, 1.)
    L            = (λ.^(1-ρ) .* (αᶠ *  Lᶠ).^ρ + (1 .- λ).^(1-ρ) .* (αᵈ * Lᵈ).^ρ).^(1/ρ)

    # Compute residuals of static equilibrium conditions
    Wages = (Wᵈ ./ Wᶠ).^(-ρ/(1-ρ)) - (Lᵈ ./ Lᶠ).^ρ .* (Domestic_Int ./ Foreign_Int).^(-ρ)
    RF    = (θ / (r + δ))^(θ / (1 - θ)) * (1 - θ) * Z .* L - (Wᵈ .* Lᵈ + Wᶠ .* Lᶠ)

    return [Wages'; RF']

end

"""
Solves the sequence of T bivariate root finding problems for each location. Returns two N × T matrices of
wages in each location.
"""
function StaticEqSolver(Transit::TransitSoln = Transit1; p::Parameters = Params1, tol::T1 = 1e-10, maxiter::T2 = Int(2e6), damper::T1 = 1e-2, inner_verbose::Bool = true) where{T1 <: Real, T2 <: Int}

    (; N) = p
    (; Wᵈ, Wᶠ, Lᵈ, Lᶠ) = Transit
    Wᵈ_new, Wᶠ_new = copy(Wᵈ), copy(Wᶠ)

    for l in 1:N - 1

        err = 1 + tol
        iter = 0
        Wᵈₗ, Wᶠₗ, Lᵈₗ, Lᶠₗ  = Wᵈ[l, :], Wᶠ[l, :], Lᵈ[l, :], Lᶠ[l, :]
        x = Wᵈₗ ./ Wᶠₗ
        while err >= tol && iter <= maxiter
            Residuals  = StaticEq(Lᵈₗ, Lᶠₗ, x, Wᶠₗ; p = p)
            err        = maximum(abs.(Residuals))
            x        .= max.(x   .+ damper .* Residuals[1, :], 1e-4)
            Wᶠₗ      .= max.(Wᶠₗ .+ damper .* Residuals[2, :], 1e-4)
            iter     += 1

            if iter % 50 == 0 && inner_verbose == true
                println("Error is $err at iter $iter")
            end
        end

        if err < tol 
            println("\tSTATIC EQUILIBRIUM: CONVERGED IN $iter ITERATIONS")
        else
            println("STATIC EQUILIBRIUM: MAX ITERATIONS REACHED AT $iter ITERATIONS")
        end

        Wᵈ_new[l, :] .= x .* Wᶠₗ
        Wᶠ_new[l, :] .= Wᶠₗ

    end

    return TransitSoln(Wᵈ_new, Wᶠ_new, Transit.Lᵈ, Transit.Lᶠ, Transit.U̇ᵈ, Transit.U̇ᶠ, Transit.Πᵈ, Transit.Πᶠ)


end

"""
This function takes the current wage vector in each location and uses it to update u̇'s
"""
function UpdateChanges(Transit::TransitSoln, μ̇::Vector{T1}; p::Parameters = Params1) where{T1 <: Real}
    
    (; νᵈ, νᶠ, N, β, Πᶠ_pre, Πᵈ_pre, T) = p
    (; U̇ᵈ, U̇ᶠ, Πᵈ, Πᶠ, Wᵈ, Wᶠ) = Transit

    # Create wage changes
    Ẇᵈ, Ẇᶠ = hcat([Wᵈ[l,t+1] / Wᵈ[l,t] for l in 1:N, t in 1:T-1], ones(N)), 
             hcat([Wᶠ[l,t+1] / Wᶠ[l,t] for l in 1:N, t in 1:T-1], ones(N))

    # Get lead of value changes - recall that boundary value has been imposed at initialization
    U̇ᵈ₊, U̇ᶠ₊ = reverse([U̇ᵈ[l,t+1] for l in 1:N, t in T-1:-1:1], dims = 2), 
               reverse([U̇ᶠ[l,t+1] for l in 1:N, t in T-1:-1:1], dims = 2)

    # Get lags of migration rates
    Πᵈ₋, Πᶠ₋ = [t == 1 ? Πᵈ_pre[l, lp] : Πᵈ[l, lp, t - 1]  for l in 1:N, lp in 1:N, t in 1:T-1],
               [t == 1 ? Πᶠ_pre[l, lp] : Πᶠ[l, lp, t - 1]  for l in 1:N, lp in 1:N, t in 1:T-1]


    # Update value changes
    U̇ᵈ_new, U̇ᶠ_new = hcat(zeros(N, T - 1), ones(N)), hcat(zeros(N, T - 1), ones(N)) # Pre-allocate, ones at the end implement the boundary for the next iteration
    for l in 1:N

        for lp in 1:N
            
            for t in 1:T - 1
                cost_change = l == N && lp < N ? μ̇[t] : 1.
                U̇ᵈ_new[l,t] += Πᵈ₋[l, lp, t] * U̇ᵈ₊[lp,t]^(β/νᵈ) * (cost_change)^(1/νᵈ)
                U̇ᶠ_new[l,t] += Πᶠ₋[l, lp, t] * U̇ᶠ₊[lp,t]^(β/νᶠ) * (cost_change)^(1/νᶠ)
            end

        end

    end

    U̇ᵈ_new = Ẇᵈ .* U̇ᵈ_new.^νᵈ
    U̇ᶠ_new = Ẇᶠ .* U̇ᶠ_new.^νᶠ

    return TransitSoln(Transit.Wᵈ, Transit.Wᶠ, Transit.Lᵈ, Transit.Lᶠ, U̇ᵈ_new, U̇ᶠ_new, Transit.Πᵈ, Transit.Πᶠ)

end

"""
Update the choice probabilities a u̇ sequence
"""
function UpdateProbabilities(Transit::TransitSoln, μ̇::Vector{T1}; p::Parameters = Params1) where{T1 <: Real}

    (; Πᵈ_pre, Πᶠ_pre, β, νᵈ, νᶠ, N, T) = p
    (; U̇ᵈ, U̇ᶠ, Πᵈ, Πᶠ) = Transit

    Πᵈ_new, Πᶠ_new = copy(Πᵈ), copy(Πᶠ)

    # Get lead of value changes - recall that boundary value has been imposed at initialization
    U̇ᵈ₊, U̇ᶠ₊ = reverse([U̇ᵈ[l,t+1] for l in 1:N, t in T-1:-1:1], dims = 2), 
               reverse([U̇ᶠ[l,t+1] for l in 1:N, t in T-1:-1:1], dims = 2)

    # Get lags of migration rates
    Πᵈ₋, Πᶠ₋ = [t == 1 ? Πᵈ_pre[l, lp] : Πᵈ[l, lp, t - 1]  for l in 1:N, lp in 1:N, t in 1:T-1],
               [t == 1 ? Πᶠ_pre[l, lp] : Πᶠ[l, lp, t - 1]  for l in 1:N, lp in 1:N, t in 1:T-1]

    # Calculate denominator and migration rates from each origin
    Denomᵈ, Denomᶠ = zeros(N, T - 1), zeros(N, T - 1) # Pre-allocate
    for l in 1:N

        for lp in 1:N
            
            for t in 1:T - 1
                cost_change = l == N && lp < N ? μ̇[t] : 1.
                Denomᵈ[l,t] += Πᵈ₋[l, lp, t] * U̇ᵈ₊[lp,t]^(β/νᵈ) * (cost_change)^(1/νᵈ)
                Denomᶠ[l,t] += Πᶠ₋[l, lp, t] * U̇ᶠ₊[lp,t]^(β/νᶠ) * (cost_change)^(1/νᶠ)
            end

        end

    end

    # Update migration rates
    for l in 1:N

        for lp in 1:N

            for t in 1:T-1
                cost_change = l == N && lp < N ? μ̇[t] : 1.
                Πᵈ_new[l,lp,t] = (U̇ᵈ₊[lp,t]^(β/νᵈ) * cost_change^(1/νᵈ) / Denomᵈ[l,t]) * Πᵈ₋[l,lp,t]
                Πᶠ_new[l,lp,t] = (U̇ᶠ₊[lp,t]^(β/νᶠ) * cost_change^(1/νᵈ) / Denomᶠ[l,t]) * Πᶠ₋[l,lp,t]
            end
        end
    end

    return TransitSoln(Transit.Wᵈ, Transit.Wᶠ, Transit.Lᵈ, Transit.Lᶠ, Transit.U̇ᵈ, Transit.U̇ᶠ, Πᵈ_new, Πᶠ_new)

end

"""
Update labor supplies.
"""
function UpdateSupplies(Transit::TransitSoln; p::Parameters = Params1)

    (; Lᵈ, Lᶠ, Πᵈ, Πᶠ) = Transit
    (; M, N, T) = p
    Lᵈ_new, Lᶠ_new = zeros(N, T), zeros(N, T)
    Lᵈ_new[:,1], Lᶠ_new[:,1] = Lᵈ[:,1], Lᶠ[:, 1]

    for t in 1:T-1

        for l in 1:N
            
            for lp in 1:N

                Lᵈ_new[l, t + 1] += Πᵈ[lp, l, t] * Lᵈ[lp, t]
                Lᶠ_new[l, t + 1] += min(Πᶠ[lp, l, t] * Lᶠ[lp, t], M[l, lp])

            end
            
        end

    end

    return TransitSoln(Transit.Wᵈ, Transit.Wᶠ, Lᵈ_new, Lᶠ_new, Transit.U̇ᵈ, Transit.U̇ᶠ, Transit.Πᵈ, Transit.Πᶠ)

end

function SolveTransition(μ̇::Vector{T1}; p::Parameters = Params1, outer_tol::T1 = 1e-6, outer_maxiter::T2 = 1000) where{T1 <: Real, T2 <: Int}

    outer_err, outer_iter = 1 + outer_tol, 0
    Soln = TransitSoln(p);

    while outer_err >= outer_tol && outer_iter <= outer_maxiter

        # Update wages
        Soln = StaticEqSolver(Soln; p = p)

        # Use the new wages to update U̇ and calculate the error and iteration counter 
        U̇ᵈ_curr, U̇ᶠ_curr = Soln.U̇ᵈ, Soln.U̇ᶠ
        Soln = UpdateChanges(Soln, μ̇; p = p)
        outer_err = max(maximum(abs.(U̇ᵈ_curr - Soln.U̇ᵈ)), maximum(abs.(U̇ᶠ_curr - Soln.U̇ᶠ)))
        outer_iter += 1

        println("****************** ITERATION $outer_iter COMLETE: Outer err = $outer_err *************************")

        # Update choice probs with the new wages
        Soln = UpdateProbabilities(Soln, μ̇ ; p = p)

        # Update labor supplies with the new choice probabilities
        Soln = UpdateSupplies(Soln; p = p)
        
    end

    return Soln
end