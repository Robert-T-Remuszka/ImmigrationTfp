import pandas as pd

Paths = {'data':'../data'}

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
                countyfip = lambda x: x['countyfip'].astype(str))
        .assign(serial = lambda x: x['serial'].str.zfill(7), # Pillow with zeros
                pernum = lambda x: x['pernum'].str.zfill(20),
                statefip = lambda x: x['statefip'].str.zfill(2),
                countyfip = lambda x: x['countyfip'].str.zfill(3))
        .assign(year = lambda x: pd.to_datetime(x['year'], 
                        format = '%Y', errors='coerce'))
        .assign(fipcode = lambda x: x['statefip'] + x['countyfip'])
        .drop(columns = ['cbserial', 'cluster','strata', 'statefip', 'countyfip'])
        
    )
    
    # Reorder columns: the combination of sample, serial and pernum uniquely identifies individuals in Ipums
    Tidy = Tidy.reindex(columns = ['sample', 'serial', 'pernum'] + 
            [col for col in Tidy.columns if col not in ['sample', 'serial', 'pernum', 'fipcode']])
    
                        
    return Tidy