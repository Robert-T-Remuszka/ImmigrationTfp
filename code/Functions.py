import numpy as np
import pandas as pd
from scipy.optimize import least_squares

Paths = {'data':'../data', 
         'cps':'../data/cps',
         'acs':'../data/acs',
         'gdp':'../data/gdp_state_industry_BEA',
         'preperiod':'../data/acs-pre-period'}

'''#################################
FUNCTIONS
#################################'''
def IpumsTidy(df):
    
    '''
    The input to this function is an ipums extract. It returns a nice and tidy dataframe
    in turn. See below for the cleaning steps.
    '''
    Tidy = (df
        .rename(columns=lambda x: x.lower())                 # Rename all columns to their lowercase counterpart
        .assign(sample = lambda x: x['sample'].astype(str),  # Change to string types
                serial = lambda x: x['serial'].astype(str),
                pernum = lambda x: x['pernum'].astype(str),
                statefip = lambda x: x['statefip'].astype(str),
                countyfip = lambda x: x['countyfip'].astype(str),
                met2013 = lambda x : x['met2013'].astype(str))
        .assign(statefip = lambda x: x['statefip'].str.zfill(2), # Pillow with zeros
                countyfip = lambda x: x['countyfip'].str.zfill(3),
                met2013 = lambda x: x['met2013'].str.zfill(5))
        .assign(year = lambda x: pd.to_datetime(x['year'], 
                        format = '%Y', errors='coerce'))
        .assign(fipcode = lambda x: x['statefip'] + x['countyfip'])
        .assign(foreign01 = lambda x: np.isin(x['citizen'], [2,3,4,5]).astype(int),  # Foreign if naturalized, not a citizen, not a citizen but has received papers, foreign born but citizen not reported
                foreign02 = lambda x: np.isin(x['citizen'], [3,4,5]).astype(int),    # Include the naturalized in the domestic workers 
                ones = 1,
                fulltime01 = lambda x: (x['uhrswork'] >= 40).astype(int),
                fulltime02 = lambda x: (x['uhrswork'] >= 35).astype(int))
    )
    
    # Reorder columns: the combination of sample, serial and pernum uniquely identifies individuals in Ipums
    Tidy = Tidy.reindex(columns = ['sample', 'serial', 'pernum', 'fipcode','met2013'] + 
            [col for col in Tidy.columns if col not in ['sample', 'serial', 'pernum', 'fipcode','met2013']])
                        
    return Tidy

def WeightedSum(df,gvars):
    
    '''
    Calculates a weighted sum of var by groups in gvars
    using weights wt. Note that the variable indicating fulltime status
    must come second in the gvars list.
    '''
    Collapsed = (
        df.groupby(gvars)
        .apply(lambda x: pd.Series({
        'HoursSupplied': np.dot(x['uhrswork'],x['perwt']),
        'BodiesSupplied': np.dot(x['ones'],x['perwt']),
        'met2013': x['met2013'].iloc[0]
        }))
        .reset_index()
        .assign(ForeignVar = gvars[0][-2:],       # Variable definition categories
                HoursVar = gvars[1][-2:])
        .rename(columns = {gvars[0]:'foreign', gvars[1]:'fulltime'})
        )
    
    return Collapsed

def BeaTidy(df):
    
    Tidy = (
        df
        .drop(columns = ['Unnamed: 0', 'CL_UNIT', 'NoteRef', 'UNIT_MULT', 'Code'])
        .rename(columns = {'GeoName':'County Name', 'DataValue':'YNom'})
        .assign(year = lambda x: pd.to_datetime(x['TimePeriod'], format='%Y', errors='coerce')) # Create a proper date var
        .drop(columns = ['TimePeriod'])                                # Drop variable
        .assign(fipcode = lambda x: x['GeoFips'].astype(str).str.zfill(5)) # Convert to string and pillow with zeros
        .drop(columns = ['GeoFips'])
    )

    return Tidy

'''#################################
OBJECTS
#################################'''
class TfpModel:
    '''
    The TfpModel object includes methods for estimating Tfp. In the current version the only estimation
    procedure implementable is Nonlinear least squares. In the future GMM may be included.

    Data must be formatted with columns [Y K F D]. See slides
    '''
    def __init__(self,Data,T,S):
        self.Data = Data                         # [logY logK F D]
        self.T = T                               # Number of years in panel
        self.S = S                               # Number of states in panel
        self.N = T * S                           # Observations
    
    def ComputeRes(self, p):
        '''
        Compute residual vector at parameter vector p. p is understood to be formatted as
        p = [S - 1 state FEs, T - 1 time FEs, task shares for each T, 4 CES parameters, Intercept]
        '''
        ρ = p[-1]  # I store the CES parameters at the end of parameter vec
        αD = p[-2] # Absolute advantage domestic
        αF = p[-3] # Absolute advantage foreign
        θ = p[-4]  # Cobb-Douglas Revenue Share
        β = p[-5]  # Intercept
        δ = np.append(0,p[:self.S - 1])                          # State fixed effects - extra 0 for conformability
        γ = np.append(0,p[self.S - 1:self.S - 1 + self.T - 1])   # Time fixed effects - extra 0 for conformability
        λ = p[self.S - 1 + self.T - 1:-5]                        # Task shares; these can be time varying too
        F, D = np.array(self.Data[:,2]), np.array(self.Data[:,3])
        
        SelectorMat = np.tile(np.identity(len(λ)), (self.S,1))
        if len(λ) == 1:
            logL = ρ **-1 * np.log((λ**(1-ρ)) * (αF * F)**ρ + ((1-λ)**(1-ρ)) * (αD * D)**ρ)
        else:
            logL = ρ**-1 * np.log((SelectorMat @ λ**(1-ρ)) * (αF * F)**ρ + (SelectorMat @ (1-λ)**(1-ρ)) * (αD * D)**ρ)
        
        logY = self.Data[:,0]
        logK = self.Data[:,1]

        # Create time fe and state fe vector
        TfeMat = np.tile(np.identity(self.T), (self.S,1))
        SfeMat = np.vstack(
            [ # Stack state indicaors in one N x (S-1) matrix; note that state S is the reference state
            np.hstack([np.ones((self.T, 1)) if i == s else np.zeros([self.T ,1]) for i in range(self.S)]) 
            for s in range(self.S)
            ]
            )
        
        return logY - (np.ones(self.N) * β + SfeMat @ δ + TfeMat @ γ + θ * logK + (1-θ) * logL)
    
    def LsEstimates(self, p0 = None):
        '''
        Estimate Tfp using nonlinear least squares
        '''
        if p0 is  None:
            Lowerbounds = (
            [-np.inf for s in range(self.S-1)] + # State fixed effects
            [-np.inf for t in range(self.T-1)] + # Time fixed effects
            [-np.inf for t in range(1)]        + # Task shares
            [-np.inf]                          + # Intercept 
            [0]                                + # CD share 
            [-np.inf for i in range(2)]        + # foreign abs advantage, domestic abs advantage
            [-np.inf])                           # CES Param
        
            Upperbounds = (
                [np.inf for s in range(self.S-1)] + # State fixed effects
                [np.inf for t in range(self.T-1)] + # Time fixed effects
                [np.inf for t in range(1)]        + # Task shares
                [np.inf]                          + # Intercept 
                [1]                               + # CD share 
                [np.inf for i in range(2)]        + # foreign abs advantage, domestic abs advantage
                [1])                                # CES Param
            
            return least_squares(self.ComputeRes, 0.5 * np.ones(self.S - 1 + self.T - 1 + 1 + 5),
                                 bounds=(Lowerbounds,Upperbounds))
        else:
            Lowerbounds = (
            [-np.inf for s in range(self.S-1)] + # State fixed effects
            [-np.inf for t in range(self.T-1)] + # Time fixed effects
            [-np.inf for t in range(self.T)]   + # Task shares
            [-np.inf]                          + # Intercept 
            [0]                                + # CD share 
            [-np.inf for i in range(2)]        + # foreign abs advantage, domestic abs advantage
            [-np.inf])                           # CES Param
        
            Upperbounds = (
                [np.inf for s in range(self.S-1)] + # State fixed effects
                [np.inf for t in range(self.T-1)] + # Time fixed effects
                [np.inf for t in range(self.T)]   + # Task shares
                [np.inf]                          + # Intercept 
                [1]                               + # CD share 
                [np.inf for i in range(2)]        + # foreign abs advantage, domestic abs advantage
                [1])                                # CES Param
        
            return least_squares(self.ComputeRes,p0,
                                 bounds=(Lowerbounds,Upperbounds))
    
    