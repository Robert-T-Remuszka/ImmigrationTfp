# The Macroeconomic Effects of Immigration: The Role of Task Specialization

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

# Run Order
There are several files that combine these raw data sources to create a panel of US states. Here are links to the files in order of which they are run and the tasks they complete;

**Remark on Raw Data:** It is not advised that you run the raw data extract codes above since all the extract output is already included in the shared data file. The extract codes are only included so that the user can see how these extracts were generated. If you would like to execute the extract codes, you will need to create a python script called ```Credentials.py``` and create a dictionary consistent with the key references in the raw download data. To do that, you will need your own API keys to the referenced APIs above. If, for some reason you find yourself running the extract code more than once, be sure to remove the previously extracted files from the location where they were saved.

1. [**Clean the Pre-Period Data**](code/CleanPrePeriod.do)
    * Output: ```data/PrePeriod.dta```
2. **API Extractions and Saving**:
    * *Remark:* The extract should be run before the read files
    * [ACS extract here](code/AcsPull.ipynb), [CPS extract here](code/CpsPull.ipynb)
    * [Read and save ACS extract](code/AcsRead.ipynb), [Read and save CPS extract](code/CpsRead.ipynb)
3. [**Clean ACS, CPS, GDP by State and Merge**](code/MakeStateAnalysisPreTfp.do)
    * Output(s): ```data/StateAnalysisPreTfp.dta```
4. [**Estimate Production Function**](code/ProdFunc_Estimate.jl)
    * [Asscoiated Types and Functions](code/ProdFunc.jl)
    * Input(s): ```data/StateAnalysisPreTfp.dta```
    * Output(s): ```data/StateTfpAndTaskAgg.csv```, ```ProductionFunction.jld2```
5. [**Merge in Production Function Outputs**](code/MakeStateAnalysis.do)
    * Input(s): ```data/StateTfpAndTaskAgg.csv```
    * Output(s): ```data/StateAnalysis.dta```
6. [**Estimate Empirical IRFs**](code/MakeIRF.do)
    * Input(s): ```data/StateAnalysis.dta```
    * Output(s): Some plots and ```IRFEstimates.dta```
7. [**Construct Initial Migration Flows**](code/MakePrePi.do)
    * Output: ```data/PiMat.dta```
8. [**Construct ACS Migration-Flow Panel**](code/MakeAcsPi.do)
    * Output: ```data/AcsPiPanel.dta```
9. [**Estimate Scale Parameters.**](code/EstimateScaleBetaAcs.do)
    * Input(s): ```data/AcsPiPanel.dta```
    * Output(s): ```data/NuBetaEstimatesAcs.dta```