import pandas as pd
import numpy as np

Paths = {'data':'../data', 
         'cps':'../data/cps',
         'acs':'../data/acs'}

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