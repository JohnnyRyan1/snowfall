#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

DESCRIPTION

See if there is a negative feedback between air temperature and snow depths.

"""

# Import modules
import numpy as np
import pandas as pd
import xarray as xr
import netCDF4
import matplotlib.pyplot as plt
from scipy import stats
import statsmodels.api as sm
from matplotlib.offsetbox import AnchoredText
import statsmodels.stats.api as sms

#%%

# Define user
user = 'johnnyryan'

# Define path
path = '/Users/' + user + '/Dropbox (University of Oregon)/research/snowfall/data/'

# Define save path
savepath = '/Users/' + user + '/Dropbox (University of Oregon)/research/snowfall/manuscript/figures/'

# Define years and regions
years = np.arange(2001, 2022)
regions = np.arange(1, 9, 1)

# Read stats DataFrame
all_stats = pd.read_csv(path + 'all_stats.csv')

#%%

for year in years:
    
    print('Processing... %.0f' %year)
    
    # Import MERRA-2 grids
    merra_sd = xr.open_dataset(path + 'merra_snowdepth/sd_' + str(year) + '.nc')
    merra_t2m = xr.open_dataset(path + 'merra_t2m/t2m_'  + str(year) + '.nc')
    
    # Find MERRA-2 grid cells that overlap with ISMIP grid
    elevation, region = [], []
    spr_temp = []
    max_snow, max_snow_day = [], []
    grid_cell_i, grid_cell_j, lon, lat = [], [], [], []
    
    # Loop over every MERRA grid cell
    for i in all_stats['grid_cell_i']:
        for j in all_stats['grid_cell_j']:
            
            # Mean Spring air temperature
            pdd = merra_t2m['t2m'][:, j, i] - 273
            pdd[pdd < 0] = 0
            temp_mean.append(np.mean(pdd[92:122]).values)
            
            # Max snow depth, max snow depth day
            max_snow.append(merra_sd['t2m'][:, j, i].values[merra_sd['t2m'][:, j, i].argmax().values])
            max_snow_day.append(merra_sd['time'][[merra_sd['t2m'][:, j, i].argmax().values]].values[0])
            
            # Grid cell
            grid_cell_i.append(i)
            grid_cell_j.append(j)
            lon.append(np.mean((min_x, max_x)))
            lat.append(np.mean((min_y, max_y)))



    # Put in DataFrame
    df = pd.DataFrame(list(zip(lon, lat, first_55_count, first_55_median, first_55_mean,
                               albedo_median, albedo_mean, elevation, region, 
                               max_snow, max_snow_day, temp_mean, 
                               grid_cell_i, grid_cell_j)))
            
    df.columns = ['lon', 'lat', 'first_55_count', 'first_55_median', 'first_55_mean', 
                  'albedo_median', 'albedo_mean', 'elevation', 'region', 
                  'max_snow', 'max_snow_day', 'temp_mean',
                  'grid_cell_i', 'grid_cell_j']
    
    # Re-index
    df.set_index(['grid_cell_i', 'grid_cell_j'], inplace=True)
    
    # Save as csv
    df.to_csv(path + 'results/results_' + str(year) + '.csv')
