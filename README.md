# Immigration, Task Specialization and Total Factor Productivity
This paper studies the timing and size of migration's productivity effects through the lens of a task-based framework. In the task based model, TFP may rise or fall in response to an exogenous flow of migration. This result has in it the reconciliation of seemingly contradictory evidence on migration's productivity effects found in the literature. I next characterize optimal migration policy in the framework and build a sufficient statistic approach from which US migration policy can be appraised. An optimal domestic migration policy sets the elasticity of output equal to the elasticity of the foreign-born consumption share of output. Implementing the framework requires estimation of the elasticities of total factor productivity, the foreign born wage, and the task aggregate with respect to an exogenous flow of migration. I build an empirical framework using the Local Projections Instrumental Variable estimator to estimate these elasticities at each horizon from one to ten years forward. Lastly, I expand the basic framework into a fully specified dynamic general equilibrium model and assess the performance of the sufficient statistic approach, using the estimated impulse response functions as target moments in a Simulated Method of Moments procedure. I use the model to assess the effects of variation in migration quotas by skill groups and solve for the optimal quotas.

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
    * Output: ```data/StateAnalysisPreTfp.dta```.
4. [**Estimate TFP**](code/EstimateTfp.jl)
    * Output: ```data/StateTfpAndTaskAgg.csv```
5. [**Merge in TFP**](code/MakeStateAnalysis.do)
    * Output: ```data/StateAnalysis.dta```
6. [**Estimate Productivity Elasticity**](code/TfpRegressions.do)
    * Output(s): ```output/graphs/IvOlsTfp.pdf```, ```output/graphs/IvLooOlsTfp.pdf```, ```FirstStageF.pdf```
