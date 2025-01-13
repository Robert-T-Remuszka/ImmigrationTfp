# Project Description

# Raw Data Sources
The data come from several sources. In order to download the data and replicate the analysis you will need your own API keys. The sources I pull from and the code that generates the raw data are:
1. [FRED](https://fred.stlouisfed.org/)
    * [Notebook File](code/FredPull.ipynb)
    * API Key Needed: Yes
2. [IPUMS USA](https://usa.ipums.org/usa/)
    * [Notebook File](code/IpumsPull.ipynb)
    * API Key Needed: Yes
3. [BEA](https://www.bea.gov/)
    * [Notebook File](code/BeaPull.ipynb)
    * API Key Needed: Yes
4. GDP by State
    * Manual download from [REAP](https://united-states.reaproject.org/)
5. Intercensal population estimates by state. Manual download from:
    1. [1990 - 1999](https://www.census.gov/data/datasets/time-series/demo/popest/intercensal-1990-2000-state-and-county-characteristics.html)
    2. [2000-2010](https://www.census.gov/data/datasets/time-series/demo/popest/intercensal-2000-2010-state.html)
    - Note: I download the "Intercensal Estiamtes of the Resident Population by Hispanic Origin"
    table since all the other tables split the population counts in more ways than this.
    3. [2010-2020](https://www.census.gov/data/datasets/time-series/demo/popest/intercensal-2010-2020-state.html)
    4. [2020-2024](https://www.census.gov/data/datasets/time-series/demo/popest/2020s-state-total.html)
    - Note: I download the table "Annual Estimates of the Resident Population for the United States, Regions
    States, District of Columbia and Puerto Rico"


# Cleaning
The code that cleans the above API pulls can be found in [MakeCountyAnalysis.ipynb](code/MakeCountyAnalysis.ipynb)

# Estimating TFP by US State
The method used for this is described in the slides. The code that implements that method can be found in
[EstimateTfp.do](code/EstimateTfp.do)

**Note:** [EstimateTfp.do](code/EstimateTfp.do) requires the command "xframeappend". If you have not
done so already then you should install this from ssc using
```Stata
ssc install xframeappend
```
