# %%####################################################################################
import pandas as pd
import numpy as np
from ipumspy import IpumsApiClient, MicrodataExtract
from Credentials import MyCredentials
from Functions import Paths as P
from pathlib import Path

CpsDir = Path(P['cps'])

# Connect to API
ipums = IpumsApiClient(MyCredentials['IpumsApiKey'])

# %% ######################################################################## PULL ASEC 1994-1999
# A list of Cps-Asec sample IDs we extract
SampleList = ['cps' + str(year) + '_03s' for year in range(1994,2000)] + ['cps2023_03s', 'cps2024_03s']

# Varibles to be included
Vars = ['STATEFIP',                              # Geographic
        'CITIZEN', 'BPL',                        # Citizenship/Nativity
        'AGE',                                   # Demographics
        'OCC', 'IND1990', 'OCC1990', 'OCC2010',  # Work
        'UHRSWORKT',                             # Work
        'EDUC',                                  # Education 
        'INCWAGE'                                # Income
        ]

for samp in SampleList:
    print('Creating and Downloading Extract: ' + samp)
    extract = MicrodataExtract('cps', [samp], Vars)
    ipums.submit_extract(extract)
    ipums.wait_for_extract(extract)
    ipums.download_extract(extract,download_dir=CpsDir)