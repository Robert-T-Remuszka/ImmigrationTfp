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
    Πᵈ₋::Matrix{T1}                             # Pre-period choice probabilities, domestic (rows origin, cols desintation)
    Πᶠ₋::Matrix{T1}                             # Pre-period choice probabilities, foreign
    Lᵈ₀::Vector{T1}                             # Time zero domestic labor supplies
    Lᶠ₀::Vector{T1}                             # Time sero foreign labor supplies

end

struct TransitSoln{T1 <: Real, T2 <: Integer}

    Wᵈ::Matrix{T1}                              # Domestic wages, rows are location and t is time
    Wᶠ::Matrix{T1}                              # Foreign wages
    Lᵈ::Matrix{T1}                              # Domestic labor supplies
    Lᶠ::Matrix{T1}                              # Foreign labor supplies
    U̇ᵈ::Matrix{T1}                              # Domestic value differences
    U̇ᶠ::Matrix{T1}                              # Foreign value changes
    Πᵈ::Array{T1}                               # Bilateral choice probabilities, domestic (rows origin, cols destination, height period)
    Πᶠ::Array{T1}                               # Bilateral choice probabilities, foreign

    T::T2                                       # Length of transition

end

#================================================================
                        CONSTRUCTOR FUNCTIONS
================================================================#
function Parameters(p::AuxParameters = p_star;
    β::T2  = 0.96,         # Equivalent to Caliendo et al quarterly discount factor of 0.99
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
    #National_Quota::T2 = (66000 + 85000) / 1e+9, # H-2B Statutory quota + H-1B Statutory with advanced degree exemption
    wᵈ_row::T2         = ( # See IRS Statistics of Income, deflated to 2009 real values using GDP deflator
    60077  / 0.77          # 1996, thousands of 2009 USD
    ) / 1e+3,
    wᶠ_row::T2         = ( # See World Bank GNI and Population Statistics
    4085 / 0.77            # 1996, thousandths of 2009 USD
    ) / 1e+3,
    νᵈ::T2 = 4.5,
    νᶠ::T2 = 4.5,
    Init_Data::DataFrame = Init_Data
    ) where{T2 <: Real}

    # How many locations are there?
    N = nrow(Init_Data);

    # Fetch the pre-period choice probabilities from the data
    foreign_cols  = filter(c -> startswith(string(c), "pi_F_"), names(Init_Data));
    domestic_cols = filter(c -> startswith(string(c), "pi_D_"), names(Init_Data));
    Πᶠ₋, Πᵈ₋ = Matrix{T2}(Init_Data[:, foreign_cols]), Matrix{T2}(Init_Data[:, domestic_cols]);

    # Fetch the time zero migration stocks - express them in millions
    Lᵈ₀, Lᶠ₀ = Vector{T2}(Init_Data[:, :Domestic_1996] ./ 1e+6), Vector{T2}(Init_Data[:, :Foreign_1996] ./ 1e+6)

    # Solve without quotas for now
    M = fill(Inf64, N, N)

    return Parameters(β, r, δ, ρ, θ, γᶠ, γᵈ, Γ, αᶠ, αᵈ, ψ, M, νᵈ, νᶠ, wᵈ_row, wᶠ_row, N, Πᵈ₋, Πᶠ₋, Lᵈ₀, Lᶠ₀)

end

function TransitSoln(p::Parameters = Params1; T0::T1 = 25) where{T1 <: Integer}

    (; N, wᵈ_row, wᶠ_row, Lᵈ₀, Lᶠ₀) = p

    # Initialize wages
    Wᵈ, Wᶠ = ones(N,T0), ones(N, T0)
    Wᵈ[N, :] .= wᵈ_row
    Wᶠ[N, :] .= wᶠ_row

    # Initialize choice probabilities --> Only require [π₀,π₁,…, π_{T-2}] --> length T-1 vector
    Πᵈ, Πᶠ = fill(1/N, (N, N, T0 - 1)), fill(1/N, (N, N, T0 - 1))

    # Initialize supplies
    Lᵈ, Lᶠ = zeros(N, T0), zeros(N, T0)
    for t in 1:T0

        for l in 1:N

            for lp in 1:N

                Lᵈ[l,t] += t == 1 ? Lᵈ₀[l] : Πᵈ[lp,l,t - 1] * Lᵈ[lp, t - 1]
                Lᶠ[l,t] += t == 1 ? Lᶠ₀[l] : Πᶠ[lp,l,t - 1] * Lᶠ[lp, t - 1]

            end

        end

    end

    # Initialize value changes - the endoenous parts are the first T-1 columns, column T always stays at 1 to enforce the boundary condition
    U̇ᵈ, U̇ᶠ = ones(N, T0), ones(N, T0)

    return TransitSoln(Wᵈ, Wᶠ, Lᵈ, Lᶠ, U̇ᵈ, U̇ᶠ, Πᵈ, Πᶠ, T0)

end

#================================================================
                        SOLVER FUNCTIONS
================================================================#
"""
Generate the residual of the relative wage and resource feasability condition for a (location, period) pair.
The inputs are u = [log(wᵈ/wᶠ), log(wᶠ)]. A root of this function is a temporary equilibrium fot (location, period)
"""
function TempEq(lᵈ::T2, lᶠ::T2, u1::T1, u2::T1; p::Parameters) where{T1 <: Real, T2 <: Real}

    (; ρ, θ, r, δ, γᶠ, γᵈ, αᶠ, αᵈ, Γ) = p
    wᶠ = exp(u2)
    wᵈ = exp(u1) * wᶠ

    Exp1           = ((γᵈ - γᶠ) * (1 - ρ)) / ((γᵈ - γᶠ) * (1 - ρ) - ρ * γᵈ)
    Exp2           = (1 - ρ) / (1 - 2 * ρ)
    base           = (αᶠ * wᵈ) / (αᵈ * wᶠ) * exp(-Γ)
    No_Arb_Stuff_F = base^((ρ * γᶠ) / ((1 - ρ) * (γᵈ - γᶠ)))
    No_Arb_Stuff_D = base^((ρ * γᵈ) / ((1 - ρ) * (γᵈ - γᶠ)))

    Foreign_Int  = max((1 - ρ) / (γᶠ * ρ) * (No_Arb_Stuff_F * Exp1 - Exp2), 1e-4)
    Domestic_Int = max(exp((ρ / (1 - ρ)) * Γ) * ((1 - ρ) / (ρ * γᵈ)) * (exp((ρ / (1-ρ)) * γᵈ) - No_Arb_Stuff_D * Exp1), 1e-4)
    Z            = max((Domestic_Int + Foreign_Int)^((1-ρ)/ρ), 1e-4)
    λ            = clamp(Foreign_Int / (Foreign_Int + Domestic_Int), 1e-4, 1.0)
    L            = (λ^(1-ρ) * (αᶠ * lᶠ)^ρ + (1 - λ)^(1-ρ) * (αᵈ * lᵈ)^ρ)^(1/ρ)

    Wages = (wᵈ / wᶠ)^(-ρ/(1-ρ)) - (lᵈ / lᶠ)^ρ * ((αᵈ^(ρ/(1-ρ)) * Domestic_Int) / (αᶠ^(ρ/(1-ρ)) * Foreign_Int))^(-ρ)
    RF    = (θ / (r + δ))^(θ / (1 - θ)) * (1 - θ) * Z * L / (wᵈ * lᵈ + wᶠ * lᶠ) - 1

    return [Wages, RF]

end

"""
Solves the sequence of T bivariate root finding problems for each location.
"""
function TempEqSolve(Transit::TransitSoln = Transit1; p::Parameters = Params1, verbose::Bool = false)

    (; N) = p
    (; Wᵈ, Wᶠ, Lᵈ, Lᶠ, T) = Transit
    Wᵈ_new, Wᶠ_new = copy(Wᵈ), copy(Wᶠ)

    for l in 1:N - 1

        for t in 1:T
            u0  = [log(Wᵈ_new[l, t] / Wᶠ_new[l, t]), log(Wᶠ_new[l, t])]
            f(u, _) = TempEq(Lᵈ[l, t], Lᶠ[l, t], u[1], u[2]; p = p)
            sol = solve(NonlinearProblem(f, u0, nothing), NewtonRaphson(); maxiters = Int(1e6), abstol = 1e-6, reltol = 1e-6)
            sol.retcode == ReturnCode.Success || error("TempEq failed at location $l, period $t (retcode=$(sol.retcode))")
            Wᵈ_new[l, t] = exp(sol.u[1]) * exp(sol.u[2])
            Wᶠ_new[l, t] = exp(sol.u[2])
        end

        verbose && println("\tSTATIC EQ location $l: done")

    end

    return TransitSoln(Wᵈ_new, Wᶠ_new, Transit.Lᵈ, Transit.Lᶠ, Transit.U̇ᵈ, Transit.U̇ᶠ, Transit.Πᵈ, Transit.Πᶠ, Transit.T)

end

"""
This function takes the current wage vector in each location and uses it to update u̇'s.
Solves the value-change recursion exactly via a backwards from the terminal condition.
"""
function UpdateChanges(Transit::TransitSoln, μ̇::Vector{T1}; p::Parameters = Params1) where{T1 <: Real}

    (; νᵈ, νᶠ, N, β, Πᶠ₋, Πᵈ₋) = p
    (; Πᵈ, Πᶠ, Wᵈ, Wᶠ, T)      = Transit

    # Initialize with boundary condition U̇[l, T] = 1
    U̇ᵈ_new = ones(N, T)
    U̇ᶠ_new = ones(N, T)

    # Backward sweep: U̇[:,t] uses freshly-computed U̇[:,t+1] rather than stale values
    for t in T-1:-1:1

        Πᵈ_t = t == 1 ? Πᵈ₋ : Πᵈ[:, :, t - 1]
        Πᶠ_t = t == 1 ? Πᶠ₋ : Πᶠ[:, :, t - 1]

        for l in 1:N
            sum_d, sum_f = 0.0, 0.0
            for lp in 1:N
                cost_change = l == N && lp < N ? μ̇[t] : 1.
                sum_d += Πᵈ_t[l, lp] * U̇ᵈ_new[lp, t+1]^(β/νᵈ) * cost_change^(-1/νᵈ)
                sum_f += Πᶠ_t[l, lp] * U̇ᶠ_new[lp, t+1]^(β/νᶠ) * cost_change^(-1/νᶠ)
            end
            ẇᵈ = t == 1 ? 1.0 : Wᵈ[l, t] / Wᵈ[l, t-1]
            ẇᶠ = t == 1 ? 1.0 : Wᶠ[l, t] / Wᶠ[l, t-1]
            U̇ᵈ_new[l, t] = ẇᵈ * sum_d^νᵈ
            U̇ᶠ_new[l, t] = ẇᶠ * sum_f^νᶠ
        end

    end

    return TransitSoln(Transit.Wᵈ, Transit.Wᶠ, Transit.Lᵈ, Transit.Lᶠ, U̇ᵈ_new, U̇ᶠ_new, Transit.Πᵈ, Transit.Πᶠ, Transit.T)

end

"""
Update the choice probabilities.
"""
function UpdateProbabilities(Transit::TransitSoln, μ̇::Vector{T1}; p::Parameters = Params1) where{T1 <: Real}

    (; Πᵈ₋, Πᶠ₋, β, νᵈ, νᶠ, N) = p
    (; U̇ᵈ, U̇ᶠ, Πᵈ, Πᶠ, T)      = Transit

    Πᵈ_new, Πᶠ_new = copy(Πᵈ), copy(Πᶠ)

    # Get lead of value changes - recall that boundary value has been imposed at initialization
    U̇ᵈ_lead, U̇ᶠ_lead = reverse([U̇ᵈ[l,t+1] for l in 1:N, t in T-1:-1:1], dims = 2),
                       reverse([U̇ᶠ[l,t+1] for l in 1:N, t in T-1:-1:1], dims = 2)

    # Get lags of migration rates
    Πᵈ_lag, Πᶠ_lag = [t == 1 ? Πᵈ₋[l, lp] : Πᵈ[l, lp, t - 1]  for l in 1:N, lp in 1:N, t in 1:T-1],
                     [t == 1 ? Πᶠ₋[l, lp] : Πᶠ[l, lp, t - 1]  for l in 1:N, lp in 1:N, t in 1:T-1]

    # Calculate denominator and migration rates from each origin
    Denomᵈ, Denomᶠ = zeros(N, T - 1), zeros(N, T - 1) # Pre-allocate
    for l in 1:N

        for lp in 1:N

            for t in 1:T - 1

                cost_change = l == N && lp != N ? μ̇[t] : 1.
                Denomᵈ[l,t] += Πᵈ_lag[l, lp, t] * U̇ᵈ_lead[lp,t]^(β/νᵈ) * (cost_change)^(-1/νᵈ)
                Denomᶠ[l,t] += Πᶠ_lag[l, lp, t] * U̇ᶠ_lead[lp,t]^(β/νᶠ) * (cost_change)^(-1/νᶠ)

            end

        end

    end

    # Update migration rates
    for l in 1:N

        for lp in 1:N

            for t in 1:T-1

                cost_change = l == N && lp != N ? μ̇[t] : 1.
                Πᵈ_new[l,lp,t] = (U̇ᵈ_lead[lp,t]^(β/νᵈ) * cost_change^(-1/νᵈ) / Denomᵈ[l,t]) * Πᵈ_lag[l,lp,t] 
                Πᶠ_new[l,lp,t] = (U̇ᶠ_lead[lp,t]^(β/νᶠ) * cost_change^(-1/νᶠ) / Denomᶠ[l,t]) * Πᶠ_lag[l,lp,t]

            end
        end
    end

    return TransitSoln(Transit.Wᵈ, Transit.Wᶠ, Transit.Lᵈ, Transit.Lᶠ, Transit.U̇ᵈ, Transit.U̇ᶠ, Πᵈ_new, Πᶠ_new, Transit.T)

end

"""
Update labor supplies.
"""
function UpdateSupplies(Transit::TransitSoln; p::Parameters = Params1)

    (; Lᵈ, Lᶠ, Πᵈ, Πᶠ, T) = Transit
    (; M, N, Lᵈ₀, Lᶠ₀)    = p
    Lᵈ_new, Lᶠ_new = zeros(N, T), zeros(N, T)
    Lᵈ_new[:,1], Lᶠ_new[:,1] = Lᵈ₀, Lᶠ₀

    for t in 1:T-1

        for l in 1:N

            for lp in 1:N

                Lᵈ_new[l, t + 1] += Πᵈ[lp, l, t] * Lᵈ[lp, t]
                Lᶠ_new[l, t + 1] += Πᶠ[lp, l, t] * Lᶠ[lp, t]

            end

        end

    end

    return TransitSoln(Transit.Wᵈ, Transit.Wᶠ, Lᵈ_new, Lᶠ_new, Transit.U̇ᵈ, Transit.U̇ᶠ, Transit.Πᵈ, Transit.Πᶠ, Transit.T)

end

"""
Extend a converged TransitSoln from T_old to T_new periods where T_old < T_new. We can then use the extended solution
as the initial sequential equilibrium for T_new > T_old.
"""
function ExtendSoln(Soln::TransitSoln, T_new::Int)

    N, T_old = size(Soln.Wᵈ)
    @assert T_new > T_old
    dT = T_new - T_old

    Wᵈ_ext  = hcat(Soln.Wᵈ, repeat(Soln.Wᵈ[:, end:end],    1, dT)) # end:end preserves the Matrix data type (i.e. recall that Wᵈ is supposed to be N×T matrix and we want to keep it that way)
    Wᶠ_ext  = hcat(Soln.Wᶠ, repeat(Soln.Wᶠ[:, end:end],    1, dT))
    U̇ᵈ_ext  = hcat(Soln.U̇ᵈ, ones(N, dT))
    U̇ᶠ_ext  = hcat(Soln.U̇ᶠ, ones(N, dT))
    Πᵈ_ext  = cat(Soln.Πᵈ,  repeat(Soln.Πᵈ[:, :, end:end], 1, 1, dT); dims = 3)
    Πᶠ_ext  = cat(Soln.Πᶠ,  repeat(Soln.Πᶠ[:, :, end:end], 1, 1, dT); dims = 3)

    Lᵈ_ext = hcat(Soln.Lᵈ, zeros(N, dT))
    Lᶠ_ext = hcat(Soln.Lᶠ, zeros(N, dT))
    Πᵈ_term, Πᶠ_term = Soln.Πᵈ[:, :, end], Soln.Πᶠ[:, :, end]
    for t in T_old:T_new-1
        for l in 1:N
            Lᵈ_ext[l, t+1] = sum(Πᵈ_term[lp, l] * Lᵈ_ext[lp, t] for lp in 1:N)
            Lᶠ_ext[l, t+1] = sum(Πᶠ_term[lp, l] * Lᶠ_ext[lp, t] for lp in 1:N)
        end
    end

    return TransitSoln(Wᵈ_ext, Wᶠ_ext, Lᵈ_ext, Lᶠ_ext, U̇ᵈ_ext, U̇ᶠ_ext, Πᵈ_ext, Πᶠ_ext, T_new)

end


"""
Different types of proportional difference sequences for μ̇. We can think of the no-shock path as
a 'factual' path. Other possible sequences can be considered 'counterfactual' paths.

******* No Shock *****
No changes in shocks so the proporional differences will always be 1.

******* AR1 Style ******
Suppose mₜ = ψmₜ₋₁ + ϵₜ. Then, if we start at the long run mean of zero (m₋₁ = 0) and hit the economy
with a shock of size m₀ at time zero and then no more shocks after
-------------------------
t | Level | Δ
-------------------------
0 | m₀    | m₀
1 | ψm₀   | -(1 - ψ)m₀
2 | ψ²m₀  | -(1 - ψ)ψm₀
3 | ψ³m₀  | -(1 - ψ)ψ²m₀
* |   *   |     *
* |   *   |     *
t | ψᵗm₀  | -(1 - ψ)ψᵗ⁻¹m₀
--------------------------
Remark: Recall that Julia's array indexing starts at one!
"""
no_shock() = t -> 1.0
ar1_shock(; m₀::Real, p::Parameters) = t -> t == 1 ? exp(m₀) : exp(-(1 - p.ψ) * p.ψ^(t - 2) * m₀)


"""
Solve for the full sequential equilibrium under a given path for μ̇
"""
function SolveSeqEq(μ̇_path::Function; p::Parameters = Params1, outer_tol::T1 = 1e-4, outer_maxiter::T2 = 10_000, init::Union{TransitSoln, Nothing} = nothing, T_step::T2 = 10, ss_tol::T1 = 1e-4) where{T1 <: Real, T2 <: Int}

    ss_err = 1 + ss_tol
    Soln   = isnothing(init) ? TransitSoln(p) : init
    μ̇      = [μ̇_path(t) for t in 1:Soln.T - 1]

    while ss_err >= ss_tol

        outer_err, outer_iter, prev_err = 1 + outer_tol, 0, Inf
        α = 0.5

        while outer_err >= outer_tol && outer_iter <= outer_maxiter

            U̇ᵈ_prev, U̇ᶠ_prev = copy(Soln.U̇ᵈ), copy(Soln.U̇ᶠ)

            Soln = UpdateProbabilities(Soln, μ̇; p = p)
            Soln = UpdateSupplies(Soln; p = p)
            Soln = TempEqSolve(Soln; p = p, verbose = outer_iter % 50 == 0)
            Soln = UpdateChanges(Soln, μ̇; p = p)
            

            # Measure the residual then take an adaptive convex combination of the new and previous Π
            outer_err = max(maximum(abs.(U̇ᵈ_prev .- Soln.U̇ᵈ)), maximum(abs.(U̇ᶠ_prev .- Soln.U̇ᶠ)))
            α         = outer_err > prev_err ? max(α * 0.5, 0.01) : min(α * 1.01, 0.9)
            U̇ᵈ_new    = α .* Soln.U̇ᵈ .+ (1 - α) .* U̇ᵈ_prev
            U̇ᶠ_new    = α .* Soln.U̇ᶠ .+ (1 - α) .* U̇ᶠ_prev
            Soln      = TransitSoln(Soln.Wᵈ, Soln.Wᶠ, Soln.Lᵈ, Soln.Lᶠ, U̇ᵈ_new, U̇ᶠ_new, Soln.Πᵈ, Soln.Πᶠ, Soln.T)
            prev_err  = outer_err
            outer_iter += 1

            if outer_iter % 50 == 0
                println("****************** ITERATION $outer_iter COMPLETE: Outer err = $outer_err (α = $α) *************************")
            end

        end

        # Recall that the T-th element of the U̇'s is already 1, so check if we get close to one in the penultimate preiod of the transition
        ss_err = max(maximum(abs.(Soln.U̇ᵈ[:, end - 1] .- 1)),
                     maximum(abs.(Soln.U̇ᶠ[:, end - 1] .- 1)))

        if ss_err >= ss_tol
            println("Extending T from $(Soln.T) to $(Soln.T + T_step) (ss_err = $ss_err)")
            T_old = Soln.T
            Soln  = ExtendSoln(Soln, T_old + T_step)
            append!(μ̇, [μ̇_path(t) for t in T_old:T_old + T_step - 1])
        end

    end

    return Soln

end

#================================================================
                        OTHER FUNCTIONS
================================================================#
function ComputeZ(Soln::TransitSoln; p::Parameters = Params1)

    (; ρ, γᶠ, γᵈ, αᶠ, αᵈ, Γ) = p
    (; Wᶠ, Wᵈ) = Soln

    # Precompute the integral expressions
    Exp1           = ((γᵈ - γᶠ) *  (1 - ρ)) / ((γᵈ - γᶠ) * (1 - ρ) - ρ * γᵈ)                       # Eω^(ρ/(1-ρ) * γᵈ / (γᵈ - γᶠ))
    Exp2           = (1 - ρ) / (1 - 2 * ρ)                                                           # Eω^(ρ / (1-ρ))
    base           = ((αᶠ .* Wᵈ) ./ (αᵈ .* Wᶠ)) .* exp(-Γ)
    No_Arb_Stuff_F = base.^((ρ * γᶠ) / ((1 - ρ) * (γᵈ - γᶠ)))                                       # γᶠ exponent for Foreign_Int  (A.9)
    No_Arb_Stuff_D = base.^((ρ * γᵈ) / ((1 - ρ) * (γᵈ - γᶠ)))                                       # γᵈ exponent for Domestic_Int (A.10)

    # Compute parts of the production function
    Foreign_Int  = max.((1 - ρ) / (γᶠ * ρ) .* (No_Arb_Stuff_F .* Exp1 .- Exp2), 1e-4)
    Domestic_Int = max.(exp((ρ / (1 - ρ)) * Γ) * ((1 - ρ) / (ρ * γᵈ)) .* (exp((ρ / (1-ρ)) * γᵈ) .- No_Arb_Stuff_D .* Exp1), 1e-4)
    Z            = max.((Domestic_Int + Foreign_Int).^((1-ρ)/ρ), 1e-4)

    return Z

end
