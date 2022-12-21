#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" 

DESCRIPTION

1. Compile years 

2. Preliminary analysis

"""

# Import modules
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Define path
path = '/Users/jryan4/Dropbox (University of Oregon)/research/snowfall/data/'

# Define years
years = np.arange(2001, 2022)

def compile_dfs(column):
    
    dfs = []
    for year in years:
        
        # Import data
        df1 = pd.read_csv(path +'results/results_' + str(year) + '.csv', index_col=['grid_cell_i', 'grid_cell_j'])
        
        # Filter grid cells with few counts
        df1 = df1[df1['first_55_count'] >= 100] # Have a look this threshold in a bit more detail
        
        # Append
        dfs.append(df1[column])
    
    # Merge
    return pd.concat(dfs, join='inner', axis=1)

snowfall = compile_dfs('snowfall_exposure')
temp = compile_dfs('temp_exposure')
expo = compile_dfs('first_55_median')

# Compute standard deviaton 
snowfall_std = snowfall.std(axis=1) / snowfall.mean(axis=1)
temp_std = temp.std(axis=1) / temp.mean(axis=1)
expo_std = expo.std(axis=1) / expo.mean(axis=1)

# Combine DataFrames
df = pd.concat([snowfall_std, temp_std, expo_std], join='inner', axis=1)
df.columns = ['snowfall', 'temp', 'expo']

"""

Describe data 

"""

# Mean timing of glacier ice exposure for a) whole ice sheet, b) by region, c) by elevation?


# Years with earlier/later glacier ice exposure?


# Where is interannual variation in glacier ice exposure higher? Regionally? Elevationally? 


# Correlations between glacier ice exposure and a) temperature and b) snowfall
























