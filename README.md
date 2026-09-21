# The Macroeconomic Effects of Immigration: The Gains from Task Specialization

# Description
The model is a dynamic discrete-choice migration model over 52 locations (50
US states + DC, referred to collectively as US "states", plus a single
aggregate "Rest of World" location) and two nativities (domestic- and
foreign-born). Workers choose locations in a forward-looking manner to maximize
expected discounted utility, with a task-specialization production side
(state-level output combining domestic and foreign labor across a continuum
of tasks) that ties migration to state-level productivity and factor shares.

# Raw Data Sources
The data come from several sources. In order to download the data and replicate the analysis you will need your own API keys. The sources I pull from and the code that generates the raw data are:
1. [**IPUMS USA**](https://usa.ipums.org/usa/)
    * [Notebook File](code/AcsPull.ipynb)
    * API Key Needed: Yes
2. [**IPUMS CPS**](https://cps.ipums.org/cps/)
    * [Notebook File](code/CpsPull.ipynb)
    * API Key Needed: Yes
3. [**GDP by State and Industry**](https://apps.bea.gov/regional/downloadzip.htm)
    * Manual download from BEA's zip file archive.
4. [**Capital Stock by State**](https://cfds.henuecon.education/index.php/data/44-yes-capital-data)
    * Method developed by El-Shagi and Yamarik (2021).
5. [**Federal Reserve Economic Data (FRED)**](https://fred.stlouisfed.org/)
    * Series: [GDP (Implicit Price Deflator, 2017 dollars)](https://fred.stlouisfed.org/series/A191RD3A086NBEA), [Gross Private Domestic Investment: Fixed Investment (Implicit Price Deflator)](https://fred.stlouisfed.org/series/A008RD3Q086SBEA)
6. [**UN Population Data**](https://population.un.org/wpp/)
7. [**UN Migrant Stock Data**](https://www.un.org/development/desa/pd/content/international-migrant-stock)
8. **World Bank Open Data** ([GDP](https://data.worldbank.org/indicator/NY.GDP.MKTP.CD), [labor force](https://data.worldbank.org/indicator/SL.TLF.TOTL.IN))
    * Fetched live, no API key needed -- see [code/FetchRowIncConstants.py](code/FetchRowIncConstants.py)
9. **IRS Statistics of Income**, Table 2 (foreign-earned income by country/region)
    * Fetched live, no API key needed -- see [code/FetchRowIncConstants.py](code/FetchRowIncConstants.py)

**Remark on Raw Data:** It is not advised that you run the raw data extract codes above since all the extract output is already included in the shared data file. The extract codes are only included so that the user can see how these extracts were generated. If you would like to execute the extract codes, you will need to create a python script called ```Credentials.py``` and create a dictionary consistent with the key references in the raw download data. To do that, you will need your own API keys to the referenced APIs above. If, for some reason you find yourself running the extract code more than once, be sure to remove the previously extracted files from the location where they were saved.

# Run Order
The pipeline runs in four phases: (1) build the state panel and estimate the
production/task-allocation parameters directly from data, (2) build the
Rest-of-World data, (3) solve the non-stochastic steady state and
calibrate the ROW linkage parameters, (4) build tables.

## Step 1: State Panel and Directly-Estimated Parameters
1. [**Clean the Pre-Period Data**](code/CleanPrePeriod.do)
    * Output: ```data/PrePeriod.dta```
2. **API Extractions and Saving**:
    * *Remark:* The extract should be run before the read files
    * [ACS extract here](code/AcsPull.ipynb), [CPS extract here](code/CpsPull.ipynb)
    * [Read and save ACS extract](code/AcsRead.ipynb), [Read and save CPS extract](code/CpsRead.ipynb)
3. [**Clean ACS, CPS, GDP by State and Merge**](code/MakeStateAnalysisPreTfp.do)
    * Output(s): ```data/StateAnalysisPreTfp.dta```
4. [**Construct ACS Migration-Flow Panel**](code/MakeAcsPi.do)
    * Output: ```data/AcsPiPanel.dta```
5. [**Construct Individual ACS Data**](code/MakeIndividualAnalysis.do)
    * Output: ```data/IndividualCpAnalysis.dta```
6. [**Estimate Scale Parameters**](code/EstimateScaleBetaAcs.do) (ν^D, ν^F)
    * Input(s): ```data/AcsPiPanel.dta```
    * Output(s): ```data/NuBetaEstimatesAcs.dta```
7. [**Estimate CA Parameters**](code/EstimateCp.do) (μ_z, ξ_z, ξ_ω)
    * Input(s): ```data/IndividualCpAnalysis.dta```
    * Output(s): ```data/CpEstimates.dta```
8. [**Estimate EOS**](code/AggSupply_Estimate.jl) (ρ, via the factor-share condition)
    * [Associated Types and Functions](code/AggSupply_Functions.jl)
    * Input(s): ```data/CpEstimates.dta```, ```data/StateAnalysisPreTfp.dta```
    * Output(s): ```data/StateTfpAndTaskAgg.csv```, ```AggSupply.jld2```
9. [**Calibrate State Capital Shares**](code/CalibrateTheta.do) (θ_l, δ_l)
    * Input(s): ```data/CapByState/state_capital_yesdata21.dta```, ```data/StateAnalysisPreTfp.dta```
    * Output(s): ```data/ThetaDelta.dta```
10. [**Merge in Production Function Outputs**](code/MakeStateAnalysis.do)
    * Input(s): ```data/StateTfpAndTaskAgg.csv```
    * Output(s): ```data/StateAnalysis.dta```
11. [**Estimate Interior Bilateral Migration Costs**](code/EstimateBilateralCosts.do) (f^n_ll' for l, l' both US states)
    * Input(s): ```data/AcsPiPanel.dta```
    * Output(s): ```data/BilateralCosts2015.dta```, ```data/AggMigrationRateByYear.dta```

## Step 2: Rest-of-World Linkage Data
These build the empirical targets and inputs needed to calibrate the ROW
migration cost and home-bias parameters, since outward flows from
the US to ROW aren't observed in the ACS/CPS.
1. [**Fetch ROW Income Data**](code/FetchRowIncConstants.py) (World Bank GDP/labor force, IRS foreign-earned income; no API keys needed)
    * Output(s): ```data/WorldBankGdpLaborForce2015.csv```, ```data/IrsForeignEarnedIncome.csv```
2. [**Construct ROW Wage Constants**](code/MakeRowIncConstants.do) (W^D_ROW, W^F_ROW, real 2009 USD)
    * Input(s): ```data/WorldBankGdpLaborForce2015.csv```, ```data/IrsForeignEarnedIncome.csv```, ```data/GdpPriceDeflator.csv```
    * Output(s): ```data/RowIncomeConstants.dta```
3. [**Construct World Population Constants**](code/MakePopulationConstants.do) (US-born abroad, total world population, both 2015)
    * Input(s): UN International Migrant Stock and World Population Prospects files (```data/UnEstimates/```)
    * Output(s): ```data/PopulationConstants.dta```
4. [**Construct ROW-Origin ACS Flow Counts**](code/MakeRowFlowRates.do) (new foreign-born immigrants, returning domestic-born natives)
    * Input(s): ```data/acs/Acs2016.dta```
    * Output(s): ```data/RowFlowCounts.dta```

## Step 3: Steady-State Solution and ROW Calibration
1. [**Steady-State Solver**](code/SsSolve.jl) ([associated types and functions](code/SsSolve_Functions.jl))
    * Defines `Parameters`, `Solution`, and the non-stochastic steady-state solve.
    * Input(s): ```AggSupply.jld2```, ```data/ThetaDelta.dta```, ```data/CpEstimates.dta```, ```data/NuBetaEstimatesAcs.dta```, ```data/BilateralCosts2015.dta```, ```data/StateAnalysisPreTfp.dta```, ```data/PopulationConstants.dta```
2. [**Calibrate ROW Costs and Home Bias**](code/RowCostsEstimate.jl) ([associated types and functions](code/RowCostsEstimate_Functions.jl)) (f^D_ROW, B^D_ROW, f^F_ROW, B^F_ROW)
    * Jointly matches two stock-share targets (share of US-born abroad, share of foreign-born in the US labor force) and two flow-rate targets (domestic return-migration rate, foreign immigration rate) via a local nonlinear solve
    * Input(s): ```data/RowIncomeConstants.dta```, ```data/PopulationConstants.dta```, ```data/RowFlowCounts.dta```, ```data/StateAnalysisPreTfp.dta```, plus everything `SsSolve_Functions.jl` needs
    * Output(s): ```code/RowCosts.jld2```

## Step 4: Tables and Other Output
1. [**Build Parameter Table**](code/MakeTables.jl)
    * Input(s): ```data/CpEstimates.dta```, ```AggSupply.jld2```, ```data/NuBetaEstimatesAcs.dta```, ```code/RowCosts.jld2```
    * Output(s): ```output/tables/ParameterEstimates.tex```