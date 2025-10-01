# Immigration, Task Specialization and Total Factor Productivity
This paper studies the timing and size of migration's productivity effects through the lens of a task-based framework. In the task based model, TFP may rise or fall in response to an exogenous flow of migration. This result has in it the reconciliation of seemingly contradictory evidence on migration's productivity effects found in the literature. I next characterize optimal migration policy in the framework and build a sufficient statistic approach from which US migration policy can be appraised. An optimal domestic migration policy sets the elasticity of output equal to the elasticity of the foreign-born consumption share of output. Implementing the framework requires estimation of the elasticities of total factor productivity, the foreign born wage, and the task aggregate with respect to an exogenous flow of migration. I build an empirical framework using the Local Projections Instrumental Variable estimator to estimate these elasticities at each horizon from one to ten years forward. Lastly, I expand the basic framework into a fully specified dynamic general equilibrium model and assess the performance of the sufficient statistic approach, using the estimated impulse response functions as target moments in a Simulated Method of Moments procedure. I use the model to assess the effects of variation in migration quotas by skill groups and solve for the optimal quotas.

# Raw Data Sources
The data come from several sources. In order to download the data and replicate the analysis you will need your own API keys. The sources I pull from and the code that generates the raw data are:
1. [FRED](https://fred.stlouisfed.org/)
    * [Notebook File](code/FredPull.ipynb)
    * API Key Needed: Yes
2. [IPUMS USA](https://usa.ipums.org/usa/)
    * [Notebook File](code/AcsPull.ipynb)
    * API Key Needed: Yes
3. [IPUMS CPS](https://cps.ipums.org/cps/)
    * [Notebook File](code/CpsPull.ipynb)
    * API Key Needed: Yes
3. [BEA](https://www.bea.gov/)
    * [Notebook File](code/BeaPull.ipynb)
    * API Key Needed: Yes
4. GDP by State and Industry
    * Manual download from BEA's [zip file archive](https://apps.bea.gov/regional/downloadzip.htm)
5. Fixed Assets by Industry
    * Manual download from BEA's [Interactive Data](https://apps.bea.gov/iTable/?ReqID=10&step=2#eyJhcHBpZCI6MTAsInN0ZXBzIjpbMiwzXSwiZGF0YSI6W1siVGFibGVfTGlzdCIsIjEyNiJdXX0=)
    * The exact account table is "Current-Cost Net Stock of Private Fixed Assets by Industry"
        - I deflate these values using the investment deflator available in fred. See the [Fred pull notebook](code/FredPull.ipynb).
6. Employment Weighted SIC/NAICS crosswalk
    * Citation: Schaller, Z., & DeCelles, P. (2021). Weighted crosswalks for NAICS and SIC industry codes. Ann Arbor, MI: Inter-university Consortium for Political and Social Research [distributor], 44-76.
    * [ICPSR Link](https://www.openicpsr.org/openicpsr/project/145101/version/V2/view;jsessionid=BA9C29B51BC66A646EDB39977723F1DB?path=/openicpsr/145101/fcr:versions/V2/Published-Crosswalk-Files&type=folder)
<!-- 
I also download population estimates. These were originally for weighting, but I ended up weighting by employment, which can be calculated from the previous data.
5. Intercensal population estimates by state. Manual download from:
    1. [1990 - 1999](https://www.census.gov/data/datasets/time-series/demo/popest/intercensal-1990-2000-state-and-county-characteristics.html)
    2. [2000-2010](https://www.census.gov/data/datasets/time-series/demo/popest/intercensal-2000-2010-state.html)
    - Note: I download the "Intercensal Estiamtes of the Resident Population by Hispanic Origin"
    table since all the other tables split the population counts in more ways than this.
    3. [2010-2020](https://www.census.gov/data/datasets/time-series/demo/popest/intercensal-2010-2020-state.html)
    4. [2020-2024](https://www.census.gov/data/datasets/time-series/demo/popest/2020s-state-total.html)
    - Note: I download the table "Annual Estimates of the Resident Population for the United States, Regions
    States, District of Columbia and Puerto Rico"
-->
# Cleaning
There are several files that combine these raw data sources to create a panel of US states. Here are links to the files in order of which they are run and the tasks they complete;
1. [Estimate Capital Stock by State](code/CapStockByStateEstimates.do)
2. [Clean ACS, CPS, GDP by State and Merge](code/MakeStateAnalysis.ipynb)
3. [Estimate and merge in TFP](code/EstimateTfp.py)

# Analysis
These files can be run in any order. The two primary files here are
1. Regressions
    - [State Level](code/TfpRegressions.do)
2. Exploration
    - [Chloropleth Map of Tfp Estimates](code/StateTfpGraphs.do)
