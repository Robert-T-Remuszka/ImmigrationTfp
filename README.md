# Project Description

# Raw Data Sources
The data come from several sources. In order to download the data and replicate the analysis you will need your own API keys. The sources I pull from and the code that generates the raw data are:
1. [FRED](https://fred.stlouisfed.org/)
    * [Notebook File](code/FredPull.ipynb)
2. [IPUMS USA](https://usa.ipums.org/usa/)
    * [Notebook File](code/IpumsPull.ipynb)
3. [BEA](https://www.bea.gov/)
    * [Notebook File](code/BeaPull.ipynb)

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
