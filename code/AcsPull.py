# %% #################################################################################### Import required libraries and set up relative paths
import pandas as pd
import numpy as np
from ipumspy import IpumsApiClient, MicrodataExtract
from Credentials import MyCredentials
from Functions import Paths as P
from pathlib import Path

DataDir = Path(P['data'])
AcsDir = Path(P['acs'])
AcsPreDir = Path(P['preperiod'])

# Connect to API
ipums = IpumsApiClient(MyCredentials['IpumsApiKey'])

# %% #################################################################################### Download and save the data
Vars = ['STATEFIP',                                                # Geographic
        'CITIZEN', 'BPL',                                          # Citizenship/nativity
        'AGE',                                                     # Demographics
        'OCC','IND1990','INDNAICS','OCC1990', 'OCC2010', 'OCCSOC', # Work
        'UHRSWORK',                                                # Work
        'EDUC',                                                    # Education
        'INCWAGE'                                                  # Income
        ]

# Only going to 2022 because 2023 doesn't have all the variables I want
SampleList =  ['us' + str(year) + 'a' for year in range(2000,2023)] 

for samp in SampleList:
    print('Creating and Downloading Extract: ' + samp)
    extract = MicrodataExtract('usa', [samp], Vars)
    ipums.submit_extract(extract)
    ipums.wait_for_extract(extract)
    ipums.download_extract(extract,download_dir=AcsDir)

# %% #################################################################################### Download "Pre-period" vars
Vars = ['STATEFIP',        # Geographic
        'BPL',             # Citizenship/nativity
        'AGE'              # Demographics 
        ]
SampleList = ['us' + str(year) + 'a' for year in range(1920,1970,10)]
for samp in SampleList:
    print('Creating and Downloading Extract: ' + samp)
    extract = MicrodataExtract('usa', [samp], Vars)
    ipums.submit_extract(extract)
    ipums.wait_for_extract(extract)
    ipums.download_extract(extract,download_dir=AcsPreDir)