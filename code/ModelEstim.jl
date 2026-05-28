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

# Choice probability matrices
labels = Init_Data[:, :State];
clim_range = (minimum(log.(Params01.Πᶠ₋)), maximum(log.(Params01.Πᶠ₋)))
p1     = Plots.heatmap(log.(Params01.Πᵈ₋), color = :viridis, xrotation = 45, xticks = (1:Params01.N, labels),
yticks = (1:Params01.N, labels), title = L"Log of $\textbf{\Pi}^D_{1995}$", xtickfontsize = 4, ytickfontsize = 4, clims = clim_range)
p2     = Plots.heatmap(log.(Params01.Πᶠ₋), color = :viridis, xrotation = 45, xticks = (1:Params01.N, labels),
yticks = (1:Params01.N, labels), title = L"Log of $\textbf{\Pi}^F_{1995}$", xtickfontsize = 4, ytickfontsize = 4, clims = clim_range)

Plots.plot(p1, p2, layout = (2,1), size = (700, 800))

# Initial distrubution of labor supplies
#exclude = ["ROW", "District of Columbia"]
#abbrevs = [state_abbrevs[s] for s in Init_Data[:, "State"] if s ∉ exclude];

## How long does it take to run the factual and counterfactual model?
