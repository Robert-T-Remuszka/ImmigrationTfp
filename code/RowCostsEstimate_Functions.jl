using DataFrames, StatFiles, JLD2, NonlinearSolve

include("SsSolve_Functions.jl")

# %% Solver Functions =========================================================

"""
Residual of the ROW-linkage calibration: given
u = (log fᵈ_ROW, log Bᵈ_ROW, log fᶠ_ROW, log Bᶠ_ROW), solves the full
non-stochastic steady state with *symmetric* ROW costs and nativity-specific
home-bias amenities Bⁿ_ROW (paper eq. 2.1-2.2's Bⁿ_c(l), normalized to 1 for
every US state per assumption A1), then returns the gap between four
model-implied moments and their targets:
  1. share of US-born people living in ROW                    (stock)
  2. 1 - Πᵈ[ROW,ROW], the domestic return-migration rate       (flow)
  3. share of foreign-born workers in the US labor force       (stock)
  4. 1 - Πᶠ[ROW,ROW], the foreign immigration rate             (flow)
A symmetric cost alone cannot move the stock share at all -- verified
directly: model-implied share_F_us was invariant to both cost and ν across a
wide grid once Wᶠ_ROW was fixed, because a symmetric cost decays both
directions of the ROW/state transition at the same rate, leaving their ratio
governed purely by the wage-driven value gap. Bⁿ_ROW breaks that
invariance by acting on the value gap directly, restoring two genuinely
independent levers per nativity: the amenity sets the stock balance, the
cost sets the flow rate's overall level.

Since flow utility is Bⁿ_c(l)·ln(w), Bⁿ_ROW enters solve_steady_state purely
through the wage argument -- Wⁿ_ROW_eff = Wⁿ_ROW_meas^Bⁿ_ROW is mathematically
exact, so no change to the core EMAX solver is needed.

Solving in logs keeps the search on the sensible cost>0, B>0 region. A
thrown error from a bad trial point is caught and turned into ones(4) (the
outer solver treats it as "a bad point" and backs off) rather than crashing.
"""
const _home_bias_call_count = Ref(0)
const _home_bias_start_time = Ref(0.0)
reset_home_bias_progress!() = (_home_bias_call_count[] = 0; _home_bias_start_time[] = time())

function home_bias_residual(u::Vector{T}, params) where {T <: Real}

    _home_bias_call_count[] += 1
    n = _home_bias_call_count[]
    elapsed = round(time() - _home_bias_start_time[], digits = 1)

    (; fixed, Wᵈ_ROW_meas, Wᶠ_ROW_meas, total_Lᵈ, total_Lᶠ,
       target_D_abroad, target_D_return, target_F_us, target_F_immigration, W0ᵈ, W0ᶠ) = params
    fᵈ_ROW, Bᵈ_ROW, fᶠ_ROW, Bᶠ_ROW = exp.(u)
    Wᵈ_ROW = Wᵈ_ROW_meas^Bᵈ_ROW
    Wᶠ_ROW = Wᶠ_ROW_meas^Bᶠ_ROW
    pstr = "fᵈ=$(round(fᵈ_ROW, sigdigits=4)) Bᵈ=$(round(Bᵈ_ROW, sigdigits=4)) " *
           "fᶠ=$(round(fᶠ_ROW, sigdigits=4)) Bᶠ=$(round(Bᶠ_ROW, sigdigits=4))"

    p = Parameters(fixed, fᵈ_ROW, fᵈ_ROW, fᶠ_ROW, fᶠ_ROW)

    try
        ss = solve_steady_state(p; Wᵈ_ROW, Wᶠ_ROW, total_Lᵈ, total_Lᶠ, W0ᵈ, W0ᶠ)
        share_D_abroad = ss.Lᵈ[p.N] / total_Lᵈ
        return_rate_D  = leave_probability(ss.Vᵈ, p.fᵈ, p.β, p.νᵈ, p.N)
        Lᵈ_us_ss = sum(ss.Lᵈ[1:p.N - 1])
        Lᶠ_us_ss = sum(ss.Lᶠ[1:p.N - 1])
        share_F_us         = Lᶠ_us_ss / (Lᵈ_us_ss + Lᶠ_us_ss)
        immigration_rate_F = leave_probability(ss.Vᶠ, p.fᶠ, p.β, p.νᶠ, p.N)
        resid = [share_D_abroad - target_D_abroad, return_rate_D - target_D_return,
                 share_F_us - target_F_us, immigration_rate_F - target_F_immigration]
        println("call $n (t=$(elapsed)s): $pstr |resid|=$(round(sqrt(sum(abs2, resid)), sigdigits=6))")
        flush(stdout)
        return resid
    catch e
        e isa InterruptException && rethrow()
        println("call $n (t=$(elapsed)s): $pstr INFEASIBLE ($(typeof(e)))")
        flush(stdout)
        return ones(T, 4)
    end

end

"""
Calibrate fᵈ_ROW, Bᵈ_ROW, fᶠ_ROW, Bᶠ_ROW to match the two stock-share and two
flow-rate targets, via NonlinearSolve.jl with a finite-difference Jacobian.
Default u0 is seeded at values a manual diagnostic search already found
well-behaved (smooth, no floating-point floors, no cost/ν invariance), so a
local TrustRegion solve converges quickly -- no grid search needed here,
unlike the abandoned asymmetric-cost/fixed-wage approach.
"""
function calibrate_home_bias(; Wᵈ_ROW_meas::Real, Wᶠ_ROW_meas::Real, total_Lᵈ::Real, total_Lᶠ::Real,
                              target_D_abroad::Real, target_D_return::Real,
                              target_F_us::Real, target_F_immigration::Real,
                              u0::Vector{<:Real} = log.([16.0, 0.978, 18.5, 1.23]),
                              W0ᵈ::Vector{<:Real} = load_wages().Wᵈ, W0ᶠ::Vector{<:Real} = load_wages().Wᶠ,
                              abstol::Real = 1e-8, maxiters::Integer = 200)

    fixed  = load_fixed_parameters()
    params = (; fixed, Wᵈ_ROW_meas, Wᶠ_ROW_meas, total_Lᵈ, total_Lᶠ,
                target_D_abroad, target_D_return, target_F_us, target_F_immigration, W0ᵈ, W0ᶠ)

    reset_home_bias_progress!()
    prob = NonlinearProblem(home_bias_residual, Float64.(u0), params)
    sol  = solve(prob, TrustRegion(autodiff = AutoFiniteDiff()); abstol, maxiters)

    sol.retcode == ReturnCode.Success || error("Home-bias calibration failed (retcode=$(sol.retcode))")

    fᵈ_ROW, Bᵈ_ROW, fᶠ_ROW, Bᶠ_ROW = exp.(sol.u)
    return (; fᵈ_ROW, Bᵈ_ROW, fᶠ_ROW, Bᶠ_ROW)

end

# %% Helper Functions =========================================================

"""
Load Wᵈ_ROW, Wᶠ_ROW (real 2009 USD) from MakeRowIncConstants.do's saved
RowIncomeConstants.dta.
"""
function load_row_income()
    d = only(eachrow(DataFrame(load(joinpath(data, "RowIncomeConstants.dta")))))
    return (Wᵈ_ROW = d.wd_row_2009usd, Wᶠ_ROW = d.wf_row_2009usd)
end

"""
Load the ACS-based ROW-origin flow counts (2015) from MakeRowFlowRates.do's
saved RowFlowCounts.dta: new foreign-born immigrants and returning
domestic-born natives, both observed in the ACS 2016 wave.
"""
function load_row_flow_counts()
    d = only(eachrow(DataFrame(load(joinpath(data, "RowFlowCounts.dta")))))
    return (immigration_count = d.immigration_count_2015, return_migration_count = d.return_migration_count_2015)
end

"""
The four ROW-cost calibration targets, all 2015: the stock share of US-born
people living abroad, the share of foreign-born workers in the US labor
force, the domestic return-migration rate, and the foreign immigration rate.
target_F_us is a US-side ratio (Lᶠ_us/(Lᵈ_us+Lᶠ_us)), not a share of world
population -- comparing US labor-force counts to total world population was
comparing incompatible units and drove the calibration to a degenerate
corner. The return/immigration rates are ACS flow counts
(load_row_flow_counts()) divided by that nativity's ROW population (total
population minus the US share) -- the empirical counterpart of the model's
1-Π^n[ROW,ROW].
"""
function row_cost_targets()

    tp    = load_total_population()
    ls    = load_labor_supply()
    pop   = only(eachrow(DataFrame(load(joinpath(data, "PopulationConstants.dta")))))
    flows = load_row_flow_counts()

    Lᵈ_us  = sum(ls.Lᵈ)
    Lᶠ_us  = sum(ls.Lᶠ)
    Lᶠ_row = tp.total_Lᶠ - Lᶠ_us
    Lᵈ_row = pop.us_abroad_2015

    target_D_abroad       = Lᵈ_row / tp.total_Lᵈ
    target_F_us           = Lᶠ_us / (Lᵈ_us + Lᶠ_us)
    target_D_return       = flows.return_migration_count / Lᵈ_row
    target_F_immigration  = flows.immigration_count / Lᶠ_row

    return (; target_D_abroad, target_F_us, target_D_return, target_F_immigration)

end

"""
Load the calibrated fᵈ_ROW, Bᵈ_ROW, fᶠ_ROW, Bᶠ_ROW from RowCostsEstimate.jl's
saved RowCosts.jld2.
"""
load_row_costs() = load(joinpath(@__DIR__, "RowCosts.jld2"), "p_star")
