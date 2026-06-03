using JLD2, StatFiles, DataFrames, TidierData, Plots, LaTeXStrings, Statistics, NonlinearSolve, LinearAlgebra
using Shapefile, GeoMakie, CairoMakie
include("Globals.jl")
include("ProdFunc.jl");
include("PF_Transition.jl");
@load "ProductionFunction.jld2" p_star;

# Get the shape file read into a DataFrame
shp_path = joinpath(data, "cb_2023_us_state_500k", "cb_2023_us_state_500k.shp");
table    = Shapefile.Table(shp_path);
shp_df   = DataFrame(table);

# Fetch the initial data
Init_Data = leftjoin(DataFrame(load(joinpath(data, "PiMat.dta"))), 
                     select(shp_df, :NAME, :geometry), 
                     on = :State => :NAME)


# Initialize some data
Params01 = Parameters();


## PLOT INITIAL DATA

labels     = Init_Data[:, :State]
tick_pos   = collect(1:Params01.N)
tick_labs  = string.(labels)
clim_range = (minimum(log.(Params01.Πᶠ₋)), maximum(log.(Params01.Πᶠ₋)))

fig = CairoMakie.Figure(size = (1600, 1400))

# Row 1: Domestic workers — heatmap then bar chart
ax1 = CairoMakie.Axis(fig[1,2], title = L"Log of $\Pi^D_{1995}$",
                       xticks = (tick_pos, tick_labs), yticks = (tick_pos, tick_labs),
                       xticklabelrotation = π/4, xticklabelsize = 5, yticklabelsize = 5)
hm  = CairoMakie.heatmap!(ax1, tick_pos, tick_pos, log.(Params01.Πᵈ₋)',
                            colormap = :viridis, colorrange = clim_range)
CairoMakie.Colorbar(fig[1,1], hm)

ax3 = CairoMakie.Axis(fig[1,3], title = "Domestic Labor Supply 1996",
                       xlabel = "Workers",
                       yticks = (tick_pos, tick_labs), yticklabelsize = 7)
CairoMakie.barplot!(ax3, tick_pos, Float64.(Init_Data[:, :Domestic_1996]), direction = :x)

# Row 2: Foreign workers — heatmap then bar chart
ax2 = CairoMakie.Axis(fig[2,2], title = L"Log of $\Pi^F_{1995}$",
                       xticks = (tick_pos, tick_labs), yticks = (tick_pos, tick_labs),
                       xticklabelrotation = π/4, xticklabelsize = 5, yticklabelsize = 5)
CairoMakie.heatmap!(ax2, tick_pos, tick_pos, log.(Params01.Πᶠ₋)',
                     colormap = :viridis, colorrange = clim_range)

ax4 = CairoMakie.Axis(fig[2,3], title = "Foreign Labor Supply 1996",
                       xlabel = L"Workers $(\log_{10})$",
                       yticks = (tick_pos, tick_labs), yticklabelsize = 7,
                       xscale = log10)
CairoMakie.barplot!(ax4, tick_pos, max.(1.0, Float64.(Init_Data[:, :Foreign_1996])),
                    direction = :x, fillto = 1.0)

CairoMakie.ylims!(ax3, 0.5, Params01.N + 0.5)
CairoMakie.ylims!(ax4, 0.5, Params01.N + 0.5)

CairoMakie.save(joinpath(graphs, "InitialData.pdf"), fig)
fig


## How long does it take to run the factual and counterfactual model?
