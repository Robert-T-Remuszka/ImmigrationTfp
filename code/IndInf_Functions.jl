#================================================================
                    PANEL LAG/LEAD CONSTRUCTION
================================================================#
"""
Return a copy of `df` with a new column `var_Lk` (k ≥ 0, a lag) or `var_Fk` (k < 0, a lead),
equal to `var` shifted k periods within each `state`. Implemented via a join on
(state, year - k) rather than row-adjacency, so genuine calendar-year gaps are respected —
matching Stata's L(k)./F(k). time-series operators under `xtset`.
"""
function add_shift(df::DataFrame, var::Symbol, k::Integer; state::Symbol = :state, year::Symbol = :year)

    newcol  = k >= 0 ? Symbol(var, "_L", k) : Symbol(var, "_F", -k)
    shifted = rename(select(df, state, year, var), var => newcol)
    shifted[!, year] = shifted[!, year] .+ k

    return leftjoin(df, shifted, on = [state, year])

end

#================================================================
                    WEIGHTED 2SLS (HANSEN FORMULA)
================================================================#
"""
Weighted two-stage least squares via the projection-matrix formula
    β̂ = (X'PW X)⁻¹X'PW y,   PW = W(W'W)⁻¹W',
applied to the weight-scaled variables (√w⊙X, √w⊙W, √w⊙y) — the vectorized IV estimator
given in Hansen's Econometrics. `X` stacks the endogenous regressor(s) and the included
exogenous controls; `W` stacks the excluded instrument(s) and the same included exogenous
controls (which serve as their own instruments), the standard 2SLS setup.
"""
function iv_2sls(y::AbstractVector, X::AbstractMatrix, W::AbstractMatrix, w::AbstractVector)

    sw = sqrt.(w)
    Xt, Wt, yt = sw .* X, sw .* W, sw .* y

    A = Wt' * Wt
    B = Wt' * Xt
    c = Wt' * yt

    β = (B' * (A \ B)) \ (B' * (A \ c))

    return vec(β)

end

#================================================================
                    SIMULATED-PANEL CONSTRUCTION
================================================================#
"""
Convert a BKM-simulated panel (`sim`, from `simulate_panel` — log-deviations from baseline)
into a DataFrame shaped for `lpiv`: `state`, `year` (the simulated period index), the four
target "levels" `Z`, `Wage_Domestic`, `Wage_Foreign`, `L` (using the baseline-normalized-to-1
convention — an additive constant/trend cancels exactly in lpiv's long difference regardless
of what it is, so treating the log-deviation itself as log(level) is exact, not approximate),
the migration-inflow impulse `fg` (one-period log growth of the foreign labor stock), and
employment weights (Baseline's own contemporaneous Lᵈ+Lᶠ over the periods overlapping the
simulated window, matching the real spec's `[pw=emp]` rather than a frozen base-period
weight). Since `fg`'s simulated variation is entirely exogenous — it's driven only by the
known shared eₜ shock, with no confounds by construction — there is no genuine endogeneity
problem in the simulated data; call `lpiv` with `instrument = :fg` here, which collapses the
2SLS estimator to OLS exactly (see lpiv's handling of instrument == fg).
"""
function panel_from_simulation(sim::NamedTuple, Baseline::Soln; p::Parameters)

    us = 1:(p.N - 1)
    N, Tsim = size(sim.Z)

    fg = fill(NaN, N, Tsim)
    fg[:, 2:end] = sim.Lᶠ[:, 2:end] .- sim.Lᶠ[:, 1:end - 1]

    df = DataFrame(
        state         = repeat(us, outer = Tsim),
        year          = repeat(1:Tsim, inner = N),
        Z             = vec(exp.(sim.Z)),
        Wage_Domestic = vec(exp.(sim.wᴰ)),
        Wage_Foreign  = vec(exp.(sim.wᶠ)),
        L             = vec(exp.(sim.L)),
        emp           = vec(Baseline.Lᵈ[us, 1:Tsim] .+ Baseline.Lᶠ[us, 1:Tsim]),
        fg            = vec(fg),
    )
    df.fg = ifelse.(isnan.(df.fg), missing, df.fg)

    return df

end

#================================================================
                    LOCAL PROJECTION - IV (LPIV)
================================================================#
"""
Estimate the paper's LPIV impulse response of `yvar` to the migration-inflow impulse `fg`,
instrumented by the contemporaneous Bartik shift-share instrument — matching
MakeIRF.do/Functions.do's `EstimateIRF` exactly: the dependent variable at horizon h is the
long difference ln(Y[t+h]) - ln(Y[t-1]); controls are `lags` lags each of the instrument, of
Y's own one-period growth (D0Y), and of `fg`; year fixed effects and `weight` (employment)
weights are included. Returns the vector of horizon-h coefficients on `fg` (h = 0,…,H), the
model/data analogue of the paper's Figure 1.
"""
function lpiv(df::DataFrame, yvar::Symbol; state::Symbol = :state, year::Symbol = :year,
    fg::Symbol = :fg, instrument::Symbol = :Bartik_1990, weight::Symbol = :emp,
    H::Integer = 9, lags::Integer = 2)

    df = add_shift(df, yvar, 1; state, year)
    ylag1 = Symbol(yvar, "_L1")
    df[!, :D0Y] = log.(df[!, yvar]) .- log.(df[!, ylag1])

    # unique()'d since instrument == fg is a valid choice (no genuine endogeneity in
    # simulated data, see panel_from_simulation) and would otherwise try to shift/select
    # the same column twice
    shift_vars = unique(vcat([(:D0Y, k) for k in 1:lags], [(fg, k) for k in 1:lags], [(instrument, k) for k in 1:lags]))
    for (v, k) in shift_vars
        df = add_shift(df, v, k; state, year)
    end

    exog_cols = unique(vcat([Symbol(:D0Y, "_L", k) for k in 1:lags],
                             [Symbol(fg, "_L", k) for k in 1:lags],
                             [Symbol(instrument, "_L", k) for k in 1:lags]))

    years   = sort(unique(df[!, year]))
    fe_cols = [Symbol("FE_", Int(yr)) for yr in years]
    for (yr, col) in zip(years, fe_cols)
        df[!, col] = Float64.(df[!, year] .== yr)
    end

    β = fill(NaN, H + 1)

    for h in 0:H

        if h == 0
            ylead = yvar
        else
            df = add_shift(df, yvar, -h; state, year)
            ylead = Symbol(yvar, "_F", h)
        end
        df[!, :Delta_h] = log.(df[!, ylead]) .- log.(df[!, ylag1])

        core_cols = unique(vcat(:Delta_h, fg, instrument, exog_cols, weight))
        sub = df[completecases(df, core_cols), vcat(core_cols, fe_cols)]

        # drop any year dummy with no variation left in this horizon's sample (rows near
        # the panel's edges get filtered out by the lag/lead structure as h grows), matching
        # how ivreg2/Stata auto-drops collinear regressors rather than erroring
        fe_active = fe_cols[[any(!=(0.0), sub[!, c]) for c in fe_cols]]

        y      = Float64.(sub[!, :Delta_h])
        Xendog = reshape(Float64.(sub[!, fg]), :, 1)
        Zinstr = reshape(Float64.(sub[!, instrument]), :, 1)
        Xexog  = reduce(hcat, [Float64.(sub[!, c]) for c in vcat(exog_cols, fe_active)])
        w      = Float64.(sub[!, weight])

        X = hcat(Xendog, Xexog)
        W = hcat(Zinstr, Xexog)

        β[h + 1] = iv_2sls(y, X, W, w)[1]   # coefficient on fg, the first column of X

    end

    return β

end
#================================================================
        INDIRECT INFERENCE: SIMULATED MOMENTS & SMD OBJECTIVE
================================================================#
const IRF_DEPVARS = [:Z, :Wage_Domestic, :Wage_Foreign, :L]

# The model is over-identified even against just two series (4 parameters vs. 2×(H+1)
# moments), and the wage IRFs were swamping estimation with dynamics the model wasn't
# getting anywhere close to regardless of θ — so only Z and L are targeted in the objective.
# simulated_moments still computes all four (IRF_DEPVARS) so the untargeted wage IRFs remain
# available to plot alongside the targeted ones.
const TARGET_DEPVARS = [:Z, :L]

"""
Compute the S-averaged simulated LPIV coefficients for all four target variables at
structural parameters (νᴰ, νᶠ, ψ) and shock scale σ, using a FIXED set of `seeds` (common
random numbers) so the resulting objective is smooth/deterministic in θ across optimizer
iterations, rather than contaminated by fresh simulation noise at every evaluation. The
expensive part — solving Baseline and the unit-shock counterfactual — happens once per
call; the S simulate+regress draws that follow are cheap by comparison.
"""
function simulated_moments(νᴰ::Real, νᶠ::Real, ψ::Real, σ::Real; Init_Data::DataFrame,
    seeds::AbstractVector{<:Integer}, Tsim::Integer = 30, K::Integer = 60, H::Integer = 9,
    T0::Integer = 100, outer_maxiter::Integer = 1000, CF_maxiter::Integer = 2000, verbose::Bool = false)

    # solve to a fixed, modest horizon rather than the true long-run steady state: the IRF
    # only needs K periods and the shock has fully decayed well before then, so with T0 giving
    # a healthy buffer past K, the boundary-condition error at T0 (discounted backward through
    # β) is negligible over the periods we actually target — ss_tol = Inf disables
    # SolveBaseline's "extend until steady state" loop entirely, converging once at T0.
    # outer_maxiter/CF_maxiter are capped well below their (10,000) library defaults so an
    # unusual/extreme θ from the outer search fails fast rather than grinding for a long time.
    p = Parameters(; Init_Data, νᵈ = νᴰ, νᶠ = νᶠ, ψ = ψ)
    Baseline = SolveBaseline(p; T0, ss_tol = Inf, outer_maxiter, verbose)
    IRF = compute_state_irf(Baseline; p, σ = 1.0, K, CF_maxiter, verbose)

    accum = Dict(v => zeros(H + 1) for v in IRF_DEPVARS)

    for seed in seeds
        sim, _ = simulate_panel(IRF, σ, Tsim; seed)
        df = panel_from_simulation(sim, Baseline; p)
        for v in IRF_DEPVARS
            accum[v] .+= lpiv(df, v; instrument = :fg, H)
        end
    end

    S = length(seeds)
    return (; (v => accum[v] ./ S for v in IRF_DEPVARS)...)

end

"""
Equal-weighted-by-variable, inverse-variance-weighted-within-variable sum of squared
deviations between the empirical LPIV coefficients (`data_β`, `data_se` — from
BaselineIRF_*.csv) and the simulated coefficients at θ, over `TARGET_DEPVARS` (Z and L only)
— an over-identified objective (4 parameters against 2×(H+1) target moments).

Within a variable, horizons are weighted by 1/Se² (so a more precisely estimated horizon
counts for more) — but each variable's within-horizon weights are separately normalized to
sum to 1 before being added together, so Z and L contribute equally to the total regardless
of the fact that Z's LPIV coefficients happen to have much smaller standard errors than L's
(~300x smaller inverse-variance on average). Without this normalization the fully efficient
diagonal weighting would make the objective essentially "match Z only" — the standard
efficient-GMM outcome, but not what's wanted here now that both series are deliberately
targeted; this equal-weight-per-variable choice trades some classical efficiency for making
sure both moment types actually discipline the estimates (in the spirit of the
identity-weighting alternative to efficient GMM discussed in Altonji and Segal (1996)).
"""
function smd_objective(νᴰ::Real, νᶠ::Real, ψ::Real, σ::Real; Init_Data::DataFrame,
    data_β::NamedTuple, data_se::NamedTuple, seeds::AbstractVector{<:Integer},
    Tsim::Integer = 30, K::Integer = 60, H::Integer = 9)

    # some corners of the search space (extreme νᴰ/νᶠ, ψ near ±1, large σ) can make the
    # nonlinear solves fail to converge or hit a domain error — treat that as "very bad"
    # rather than letting one candidate crash the whole batch (e.g. pmap over many candidates)
    sim_β = try
        simulated_moments(νᴰ, νᶠ, ψ, σ; Init_Data, seeds, Tsim, K, H)
    catch e
        @warn "smd_objective: solve failed at θ=($νᴰ, $νᶠ, $ψ, $σ): $e — returning a large penalty"
        return 1e10
    end

    Q = 0.0
    for v in TARGET_DEPVARS
        raw_w = 1 ./ data_se[v] .^ 2
        w = raw_w ./ sum(raw_w)
        Q += sum(w .* (data_β[v] .- sim_β[v]) .^ 2)
    end

    return Q

end
