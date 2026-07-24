# %% Setup — spin up worker processes first, then load packages/includes on all of them.
# Nelder-Mead's sequential nature wastes almost all of a multi-core allocation, so we instead
# do a parallel global search (many candidate θ's evaluated simultaneously, one per worker)
# followed by a quick sequential local refinement from the best point found.
using Distributed

# prefer the SLURM allocation's own count over Sys.CPU_THREADS, since it isn't guaranteed
# that a cgroup-restricted job correctly reports the allocated (rather than the node's full)
# core count through that query
N_CORES = something(
    tryparse(Int, get(ENV, "SLURM_CPUS_PER_TASK", "")),
    tryparse(Int, get(ENV, "SLURM_JOB_CPUS_PER_NODE", "")),
    Sys.CPU_THREADS,
)
N_WORKERS = max(N_CORES - 2, 1)   # leave 2 cores for the driver process + OS overhead
addprocs(N_WORKERS)

@everywhere using JLD2, StatFiles, DataFrames, CSV, NonlinearSolve, LinearAlgebra, ForwardDiff, Random

# force single-threaded BLAS per worker — with this many worker processes, letting each
# also spawn multiple BLAS threads would oversubscribe the machine and slow things down; the
# model's linear algebra (51 locations) is far too small to benefit from BLAS threading anyway
@everywhere LinearAlgebra.BLAS.set_num_threads(1)

# absolute paths (anchored to this file's own location), since worker processes spawned by
# addprocs() are not guaranteed to share the launching process's working directory
@everywhere const CODE_DIR = @__DIR__
@everywhere include(joinpath(CODE_DIR, "Globals.jl"))
@everywhere include(joinpath(CODE_DIR, "ProdFunc.jl"))
@everywhere include(joinpath(CODE_DIR, "Solve_Baseline_Functions.jl"))
@everywhere include(joinpath(CODE_DIR, "Solve_Counterfactual_Functions.jl"))
@everywhere include(joinpath(CODE_DIR, "Simulation_Functions.jl"))
@everywhere include(joinpath(CODE_DIR, "IndInf_Functions.jl"))

using Optim, Plots

# %% Load the empirical target (coefficients + SEs) that indirect inference will match
Init_Data = load_init_data()
depvars = [:Z, :Wage_Domestic, :Wage_Foreign, :L]

data_β  = (; (v => CSV.read(joinpath(data, "BaselineIRF_$(v).csv"), DataFrame).Beta_Iv1990 for v in depvars)...)
data_se = (; (v => CSV.read(joinpath(data, "BaselineIRF_$(v).csv"), DataFrame).Se_Iv1990   for v in depvars)...)

S_estimation = 20
seeds_estimation = collect(1:S_estimation)

# νᴰ, νᶠ, ψ, σ — νᴰ/νᶠ tightened to CDP's own ballpark (≈4-5) per the user, kept away from
# the true extremes of what's structurally valid: |ψ|→1 is near-unit-root, and BKM's own
# linearization logic wants σ (the shock size) genuinely "small," so σ=5 (an e⁵≈148x
# mobility-cost swing) is well outside anything plausible and was very likely what was
# causing solves to hang rather than fail cleanly
lower = [0.5, 0.5, -0.9, 0.05]
upper = [6.0, 6.0, 0.9, 2.0]

# Broadcast the shared (read-only) data to every worker, then define the objective as an
# @everywhere function rather than a local closure — pmap can ship plain data (candidates)
# to a function each worker already has compiled locally, but Julia 1.12 will not reliably
# ship a fresh closure's method definition to workers on demand (confirmed: this errors with
# "UndefVarError: #eval_θ not defined in Main" if eval_θ is a local closure instead).
@everywhere Init_Data = $Init_Data
@everywhere data_β = $data_β
@everywhere data_se = $data_se
@everywhere seeds_estimation = $seeds_estimation

@everywhere eval_θ(θ) = smd_objective(θ[1], θ[2], θ[3], θ[4]; Init_Data, data_β, data_se, seeds = seeds_estimation)

# %% Parallel global search: draw N_global candidate θ's spanning the bounded box and
# evaluate them ALL in parallel, one per worker per "wave" — each evaluation still costs
# ~2 minutes (its own baseline + counterfactual solve), but with N_WORKERS this many, the
# whole global search finishes in only a couple of waves instead of N_global sequential steps.
N_global = 4 * N_WORKERS
Random.seed!(2026)
candidates = [lower .+ rand(4) .* (upper .- lower) for _ in 1:N_global]

println("Launching parallel global search: $N_global candidates across $N_WORKERS workers")
flush(stdout)
Qs = pmap(eval_θ, candidates)

best_idx = argmin(Qs)
θ_best = candidates[best_idx]
println("Best from global search: θ = $θ_best   Q = $(Qs[best_idx])")
flush(stdout)

# %% Local refinement from the best global point (sequential — now polishing near a good
# basin rather than exploring, so this should converge quickly).
# Abandoned Fminbox's barrier approach after it let the inner NelderMead step to ψ=-1.28
# against a [-0.9,0.9] bound (confirmed 2026-07-24) — the barrier's influence near a boundary
# depends on μ being large enough relative to the objective's own scale (Q~950), and tuning
# that reliably is more fragile than just removing the failure mode outright.
# A pure clamp-then-evaluate wrapper (tried first) is NOT enough on its own: clamping creates
# a flat plateau wherever multiple raw points map to the same clamped value, which confuses
# NelderMead's simplex geometry — verified this directly: on a toy quadratic with its true
# minimum safely *inside* the bounds, a pure clamp wrapper still converged to the boundary
# rather than the true interior optimum. Fix: clamp before the (expensive, crash-prone)
# evaluation for safety, but ALSO add a smooth quadratic exterior penalty for how far outside
# the box the raw (unclamped) point is — this gives NelderMead a real gradient-like signal to
# correct itself with, rather than a flat plateau, while still guaranteeing every actual
# smd_objective call happens at a feasible point. Verified on the same toy problem: correctly
# finds the true interior minimum when unconstrained, and correctly settles on the boundary
# when the true optimum lies outside it.
function eval_θ_bounded(θ)
    penalty = sum(max.(0.0, lower .- θ) .^ 2) + sum(max.(0.0, θ .- upper) .^ 2)
    return eval_θ(clamp.(θ, lower, upper)) + 1e4 * penalty
end

result = Optim.optimize(eval_θ_bounded, θ_best, NelderMead(),
    Optim.Options(time_limit = 3 * 3600.0, show_trace = true, show_every = 1))

θ_local = clamp.(Optim.minimizer(result), lower, upper)
Q_local = eval_θ(θ_local)

θ̂, Q̂ = Q_local < Qs[best_idx] ? (θ_local, Q_local) : (θ_best, Qs[best_idx])
νᴰ̂, νᶠ̂, ψ̂, σ̂ = θ̂
println("Estimated θ̂: νᴰ=$νᴰ̂, νᶠ=$νᶠ̂, ψ=$ψ̂, σ=$σ̂   (Q=$Q̂)")

jldsave(joinpath(data, "IndInfEstimate.jld2"); θ̂, νᴰ̂, νᶠ̂, ψ̂, σ̂, Q = Q̂)

# %% Plot the estimated model's IRF against the real Figure-1 LPIV coefficients (with the
# empirical ±1.645·SE band) — the key output of this estimation. Wrapped in the same
# try/catch pattern as smd_objective for defense in depth (this direct call bypasses that
# protection otherwise — the uncaught version of this exact call is what crashed the job).
sim_β_hat = try
    simulated_moments(νᴰ̂, νᶠ̂, ψ̂, σ̂; Init_Data, seeds = seeds_estimation)
catch e
    @warn "Final simulated_moments call failed at the estimated θ̂: $e — skipping the plot"
    nothing
end

if !isnothing(sim_β_hat)

    # Z, L (the targeted moments) share the top row; the two untargeted wage IRFs share the
    # bottom row. Model line is orange throughout — solid where targeted, dashed where not.
    plot_order = [:Z, :L, :Wage_Domestic, :Wage_Foreign]
    ylabs = Dict(:Z => "Z", :Wage_Domestic => "wᴰ", :Wage_Foreign => "wᶠ", :L => "L")
    plots = map(plot_order) do v
        h = 0:9
        lo = data_β[v] .- 1.645 .* data_se[v]
        hi = data_β[v] .+ 1.645 .* data_se[v]
        targeted = v in TARGET_DEPVARS
        plt = plot(h, data_β[v], ribbon = (data_β[v] .- lo, hi .- data_β[v]),
            fillalpha = 0.2, label = "Data", linewidth = 2, color = :steelblue,
            xlabel = "h", ylabel = "Δʰln($(ylabs[v]))", title = string(v), grid = false)
        plot!(plt, h, sim_β_hat[v], label = targeted ? "Model (targeted)" : "Model (untargeted)",
            linewidth = 2, color = :orange, linestyle = targeted ? :solid : :dash)
        plt
    end

    irf_comparison = plot(plots..., layout = (2, 2), size = (1000, 700))
    savefig(irf_comparison, joinpath(graphs, "ModelVsDataIRF.pdf"))
    irf_comparison

end