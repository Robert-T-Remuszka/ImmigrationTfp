using DataFrames, StatFiles, JLD2, NonlinearSolve

include("SsSolve_Functions.jl")

# %% Solver Functions =========================================================

"""
Residual of the ROW calibration: given
u = (log fᵈ_ROW, log wᵈ_ROW, log fᶠ_ROW, log wᶠ_ROW), solves the full
steady state with *symmetric* ROW costs, calibrated Rest of World wages wⁿ_ROW
and the entry lottery with admission probability ξ (held fixed), then returns
the gap between four model-implied moments and their targets:
  1. share of US-born people living in ROW                    (stock)
  2. 1 - Πᵈ[ROW,ROW], the domestic return-migration rate       (flow)
  3. share of foreign-born workers in the US labor force       (stock)
  4. ξ·(1 - π̃ᶠ[ROW,ROW]), the realized foreign immigration rate (flow)
For each nativity, the symmetric cost sets the flow rate's overall level and
the Rest of World wage sets the stock balance through the value gap. ξ does
not enter the application probabilities π̃, only the realized flows.

Solving in logs keeps the search on the cost>0, w>0 region. A thrown error
from a bad trial point is caught and turned into ones(4) (the outer solver
treats it as "a bad point" and backs off) rather than crashing.
"""
const _row_calibration_call_count = Ref(0)
const _row_calibration_start_time = Ref(0.0)
reset_row_calibration_progress!() = (_row_calibration_call_count[] = 0; _row_calibration_start_time[] = time())

"""
The model-implied calibration moments at a solved steady state `ss`: the four
targeted in the residual and the foreign-born application share 1 - π̃ᶠ[ROW,ROW].
"""
function row_calibration_moments(ss::Solution, p::Parameters, ξ::Real, total_Lᵈ::Real)
    L = p.N - 1
    share_D_abroad      = ss.Lᵈ[p.N] / total_Lᵈ
    return_rate_D       = leave_probability(ss.Vᵈ, p.fᵈ, p.β, p.νᵈ, p.N)
    share_F_us          = sum(ss.Lᶠ[1:L]) / (sum(ss.Lᵈ[1:L]) + sum(ss.Lᶠ[1:L]))
    application_share_F = leave_probability(ss.Vᶠ, p.fᶠ, p.β, p.νᶠ, p.N)
    immigration_rate_F  = ξ * application_share_F
    return (; share_D_abroad, return_rate_D, share_F_us, immigration_rate_F, application_share_F)
end

function row_calibration_residual(u::Vector{T}, params) where {T <: Real}

    _row_calibration_call_count[] += 1
    n = _row_calibration_call_count[]
    elapsed = round(time() - _row_calibration_start_time[], digits = 1)

    (; fixed, ξ, total_Lᵈ, total_Lᶠ,
       target_D_abroad, target_D_return, target_F_us, target_F_immigration, W0ᵈ, W0ᶠ) = params
    fᵈ_ROW, wᵈ_ROW, fᶠ_ROW, wᶠ_ROW = exp.(u)
    pstr = "fᵈ=$(round(fᵈ_ROW, sigdigits=4)) wᵈ=$(round(wᵈ_ROW, sigdigits=5)) " *
           "fᶠ=$(round(fᶠ_ROW, sigdigits=4)) wᶠ=$(round(wᶠ_ROW, sigdigits=5))"

    p = Parameters(fixed, fᵈ_ROW, fᵈ_ROW, fᶠ_ROW, fᶠ_ROW)

    try
        ss = solve_steady_state(p; Wᵈ_ROW = wᵈ_ROW, Wᶠ_ROW = wᶠ_ROW, total_Lᵈ, total_Lᶠ, W0ᵈ, W0ᶠ, ξ)
        m = row_calibration_moments(ss, p, ξ, total_Lᵈ)
        resid = [m.share_D_abroad - target_D_abroad, m.return_rate_D - target_D_return,
                 m.share_F_us - target_F_us, m.immigration_rate_F - target_F_immigration]
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
Calibrate fᵈ_ROW, wᵈ_ROW, fᶠ_ROW, wᶠ_ROW to match the two stock-share and two
flow-rate targets at a fixed admission probability ξ, via NonlinearSolve.jl
with a finite-difference Jacobian. `u0` is the log of the starting point;
for small ξ, start from the solution at a nearby larger ξ.
"""
function calibrate_row(; ξ::Real, total_Lᵈ::Real, total_Lᶠ::Real,
                        target_D_abroad::Real, target_D_return::Real,
                        target_F_us::Real, target_F_immigration::Real,
                        u0::Vector{<:Real},
                        W0ᵈ::Vector{<:Real} = load_wages().Wᵈ, W0ᶠ::Vector{<:Real} = load_wages().Wᶠ,
                        abstol::Real = 1e-8, maxiters::Integer = 200)

    fixed  = load_fixed_parameters()
    params = (; fixed, ξ = Float64(ξ), total_Lᵈ, total_Lᶠ,
                target_D_abroad, target_D_return, target_F_us, target_F_immigration, W0ᵈ, W0ᶠ)

    reset_row_calibration_progress!()
    prob = NonlinearProblem(row_calibration_residual, Float64.(u0), params)
    sol  = solve(prob, TrustRegion(autodiff = AutoFiniteDiff()); abstol, maxiters)

    sol.retcode == ReturnCode.Success || error("ROW calibration failed at ξ=$ξ (retcode=$(sol.retcode))")

    fᵈ_ROW, wᵈ_ROW, fᶠ_ROW, wᶠ_ROW = exp.(sol.u)
    return (; fᵈ_ROW, wᵈ_ROW, fᶠ_ROW, wᶠ_ROW)

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
domestic-born natives (full-time workers, 16+), and new foreign-born
immigrants aged 15+ of any work status, all observed in the ACS 2016 wave.
"""
function load_row_flow_counts()
    d = only(eachrow(DataFrame(load(joinpath(data, "RowFlowCounts.dta")))))
    return (immigration_count = d.immigration_count_2015, return_migration_count = d.return_migration_count_2015,
            adult_immigration_count = d.adult_immigration_count_2015)
end

"""
The five ROW calibration targets, all 2015: the stock share of US-born
people living abroad, the share of foreign-born workers in the US labor
force, the domestic return-migration rate, the foreign immigration rate, and
the foreign-born application share.
target_F_us is a US-side ratio (Lᶠ_us/(Lᵈ_us+Lᶠ_us)), not a share of world
population -- comparing US labor-force counts to total world population was
comparing incompatible units and drove the calibration to a degenerate
corner. The return/immigration rates are ACS flow counts
(load_row_flow_counts()) divided by that nativity's ROW population (total
population minus the US share) -- the empirical counterpart of the model's
1-Π^n[ROW,ROW]. The application share is the number of people who would like
to move to the US over the foreign-born ROW population, the counterpart of
the model's 1-π̃ᶠ[ROW,ROW]. Gallup World Poll, 2015-2017 interviews: about
158 million adults worldwide name the US as their desired destination ("More
Than 750 Million Worldwide Would Migrate If They Could", Dec 10, 2018,
https://news.gallup.com/poll/245255/750-million-worldwide-migrate.aspx).
"""
const gallup_us_desired = 158e6

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
    target_F_application  = gallup_us_desired / Lᶠ_row

    return (; target_D_abroad, target_F_us, target_D_return, target_F_immigration, target_F_application)

end

"""
Load the home-bias calibration (fᵈ_ROW, Bᵈ_ROW, fᶠ_ROW, Bᶠ_ROW; 51 US regions,
no quota) saved in RowCosts.jld2. The current model's transition numbers use it.
"""
load_row_costs() = load(joinpath(@__DIR__, "RowCosts.jld2"), "p_star")

"""
Load the ROW calibration with the entry lottery from RowCostsEstimate.jl's
saved RowCostsQuota.jld2: `which = "baseline"` (ξ from Gallup) or
`"no_quota"` (ξ = 1).
"""
load_row_costs_quota(which::AbstractString = "baseline") = load(joinpath(@__DIR__, "RowCostsQuota.jld2"), which)
