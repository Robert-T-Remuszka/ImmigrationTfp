# Immigration, Task Specialization and Total Factor Productivity
This paper studies migration's productivity effects through the lens of a task-based framework.

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

# Data Preparation
There are several files that combine these raw data sources to create a panel of US states. Here are links to the files in order of which they are run and the tasks they complete;

**Remark on Raw Data:** Extract codes are included so that you can see how these extracts were generated. If you would like to execute the extract codes, you will need to create a python script called ```Credentials.py``` and create a dictionary consistent with the key references in the raw download data. To do that, you will need your own API keys to the referenced APIs above.

1. [**Clean the Pre-Period Data**](code/CleanPrePeriod.do)
    * Output: ```data/PrePeriod.dta```
2. **API Extractions and Saving**:
    * *Remark:* The extract should be run before the read files
    * [ACS extract here](code/AcsPull.ipynb), [CPS extract here](code/CpsPull.ipynb)
    * [Read and save ACS extract](code/AcsRead.ipynb), [Read and save CPS extract](code/CpsRead.ipynb)
3. [**Clean ACS, CPS, GDP by State and Merge**](code/MakeStateAnalysisPreTfp.do)
    * Output: ```data/StateAnalysisPreTfp.dta```.
4. [**Estimate Production Function**](code/ProdFunc_Estimate.jl)
    * Internal Dependencies: [ProdFunc.jl](code/ProdFunc.jl)
    * Output: ```data/StateTfpAndTaskAgg.csv```
5. [**Merge in Production Function Output**](code/MakeStateAnalysis.do)
    * Output: ```data/StateAnalysis.dta```
6. [**Construct Migration Flows**](code/MakePrePi.do)
    * Output: ```data/PiMat.dta```

# Analysis Files

1. [**Create Empirical IRFs**](code/MakeIRF.do)
2. [**Estimate Dynamic Migration Parameters**](code/IndInf_Estimate.jl)
    * Internal Dependencies: [Solve_Baseline_Functions.jl](code/Solve_Baseline_Functions.jl),