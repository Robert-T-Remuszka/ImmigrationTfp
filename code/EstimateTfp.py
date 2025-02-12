# %%
import numpy as np
import pandas as pd
from Functions import *

# %%
StateLong = pd.read_stata(Paths['data'] + '/StateAnalysisFile.dta')
leaveout = ['ImmigrantGroup','foreign','statefip','year','HoursSupplied','BodiesSupplied','Wage','StateName',
            'NGdp','KNom','PriceDeflator','InvestmentDeflator']

StateAggregates = {var:'mean' for var in StateLong.columns if var not in leaveout} | {'Y':'mean', 'K':'mean'}

StateLong = (
    pd.read_stata(Paths['data'] + '/StateAnalysisFile.dta')
    .assign(
        Y = lambda x: (x['NGdp'] * 100 / x['PriceDeflator']) * 1e+6, # Units of gdp WERE millions of USD
        # These were likely missing CITIZEN, or born in US territories other than PR
        foreign = lambda x: x['foreign'].mask((x['foreign'] == 0) & (x['ImmigrantGroup'] != 'United States') , 1))
    .rename(columns =  {'K':'KNom'})
    .assign(K = lambda x: (x['KNom'] * 100 / x['InvestmentDeflator']))
    .groupby(['statefip','StateName','year','foreign'])
    .agg({'BodiesSupplied':'sum'} | StateAggregates)
    .reset_index()
    .assign(logY = lambda x: np.log(x['Y']),
            logK = lambda x: np.log(x['K']),
            year = lambda x: x['year'].dt.year.astype(int))

)

Foreign = (
    StateLong.loc[StateLong['foreign'] == 1]
    .rename(columns={'BodiesSupplied':'F'})
    .drop(columns = ['foreign'])
)
Domestic = (
    StateLong.loc[StateLong['foreign'] == 0]
    .rename(columns={'BodiesSupplied':'D'})
    .pipe(lambda x: x.loc[:,['statefip','StateName','year','D']])
)
AnalysisDf = Foreign.merge(Domestic,how='left',on=['statefip','StateName','year'])

# %%
S, T = 51, 2023 - 1994 + 1
Data = AnalysisDf[['logY', 'logK', 'F', 'D']].to_numpy()
'''
A model without any time variation in the task shares can be estimated using,

TfpModelObj = TfpModel(Data, T, S).LsEstimates(p0=np.hstack([np.ones(S - 1 + T - 1), # State and time Fes
                                                            0.50 * np.ones(1), # Task shares
                                                            0.85 * np.ones(3), # Intercept, absolute advantages
                                                            0.30 * np.ones(1),  # Capital share
                                                            0.85 * np.ones(1)])) # CES param

This was the version of the model used in earlier versions of the paper to generate TFP estimates.
'''
TfpModelObj = TfpModel(Data, T, S).LsEstimates(p0=np.hstack([np.ones(S - 1 + T - 1), # State and time Fes
                                                            0.50 * np.ones(T), # Task shares
                                                            np.ones(3), # Intercept, absolute advantages
                                                            0.30 * np.ones(1), # Capital share
                                                            0.85 * np.ones(1)])) # CES param


# %%
θ = TfpModelObj.x[-4]  # Cobb-Douglas Revenue Share
β = TfpModelObj.x[-5]  # Intercept
δ = np.append(0,TfpModelObj.x[:S - 1]) # State fixed effects
SfeMat = np.vstack( # A useful matrix for calculating TfpEstimates
            [ # Stack state indicaors in one N x (S-1) matrix; note that state S is the reference state
            np.hstack([np.ones((T, 1)) if i == s else np.zeros([T,1]) for i in range(S)]) 
            for s in range(S)
            ]
            )
AnalysisDf['Z'] = np.exp((TfpModel(Data, T, S).ComputeRes(TfpModelObj.x) + SfeMat @ δ + β) * (1/(1-θ)))

# %%
# Saving
AnalysisDf = AnalysisDf.rename(columns={'D':'BodiesSupplied0', 'F':'BodiesSupplied1'})
AnalysisDf.to_stata(Paths['data'] + '/StateAnalysisFileTfp.dta', write_index=False)