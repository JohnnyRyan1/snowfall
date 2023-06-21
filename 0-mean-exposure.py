#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" 

DESCRIPTION

1. Find min, mean, max glacier ice exposure date for each region and by elevation.

"""

# Import modules
import xarray as xr
import netCDF4
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Define path
path = '/Users/jryan4/Dropbox (University of Oregon)/research/snowfall/data/'
path_alt = '/Users/jryan4/Dropbox (University of Oregon)/research/clouds/data/'

# Import regions
regions_file = netCDF4.Dataset(path_alt + 'temp_albedo_summer_climatologies.nc')
regions = regions_file.variables['regions'][:].filled()
summer_temp = regions_file.variables['t2m_mean'][:].filled()

# Import snowfall climatlogy
snowfall_file = netCDF4.Dataset(path + 'cloudsat/cloudsat-greenland-snowfall-2006-2016.nc')
snowfall = snowfall_file.variables['mean_annual_snowfall'][:].filled()

# Import ISMIP 1 km grid
ismip_1km = xr.open_dataset(path + 'masks/1km-ISMIP6.nc')

# Define years
years = np.arange(2001, 2022)

regional_median_55 = []
regional_median_60 = []
regional_count_55 = []
regional_count_60 = []
regional_temp = []

for year in years:
    
    print('Processing... %.0f' %year)
       
    # Import MODIS grids
    modis = xr.open_dataset(path + 'modis_exposure_v2/MOD10A1_albedo_dates_' + str(year) + '.nc')
    
    # Produce glacier ice masks
    ice_mask_first_55 = modis['first_55'].values.astype(float)
    ice_mask_first_55[ismip_1km['ICE'] == 0] = np.nan
    ice_mask_first_55[ice_mask_first_55 == 0] = np.nan
    
    ice_mask_first_60 = modis['first_60'].values.astype(float)
    ice_mask_first_60[ismip_1km['ICE'] == 0] = np.nan
    ice_mask_first_60[ice_mask_first_60 == 0] = np.nan

    ice_mask_last_55 = modis['last_55'].values.astype(float)
    ice_mask_last_55[ismip_1km['ICE'] == 0] = np.nan
    ice_mask_last_55[ice_mask_last_55 == 0] = np.nan

    ice_mask_last_60 = modis['last_60'].values.astype(float)
    ice_mask_last_60[ismip_1km['ICE'] == 0] = np.nan
    ice_mask_last_60[ice_mask_last_60 == 0] = np.nan
    
    ice_mask_elevation = ismip_1km['SRF'].values
    ice_mask_elevation[ice_mask_elevation == 0] = np.nan

    data_55 = []
    data_60 = []
    count_55 = []
    count_60 = []
    temp = []
    
    # Group by elevation and region
    for e in np.arange(0, 3550, 100):
        for r in np.arange(1, 9, 1):
            #print(e, r)
            mask = (np.isfinite(ice_mask_first_55) & (regions == r) & (ice_mask_elevation == e))
            data_55.append(np.nanmedian(ice_mask_first_55[mask]))
            data_60.append(np.nanmedian(ice_mask_first_60[mask]))
            count_55.append(np.nansum(ice_mask_first_55[mask]))
            count_60.append(np.nansum(ice_mask_first_60[mask]))
            temp.append(np.nanmean(summer_temp[mask]))
            
    # Reshape and append
    regional_median_55.append(np.reshape(np.array(data_55), (36, 8)))
    regional_median_60.append(np.reshape(np.array(data_60), (36, 8)))
    regional_count_55.append(np.reshape(np.array(count_55), (36, 8)))
    regional_count_60.append(np.reshape(np.array(count_60), (36, 8)))
    regional_temp.append(np.reshape(np.array(temp), (36, 8)))
                              
# Median glacier ice exposure by region (independent of elevation)
df_median = pd.DataFrame(list(zip(np.nanmedian(np.nanmedian(np.array(regional_median_55), axis=1), axis=0))),
                         ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW'])
df_median.columns = ['exposure']

# Mean glacier ice exposure by elevation for each region?
df = pd.DataFrame((np.nanmean(np.array(regional_median_55)[:,:,:], axis=0)))
df.columns = ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW']

df_count = pd.DataFrame((np.nanmean(np.array(regional_count_55)[:,:,:], axis=0)))
df_count.columns = ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW']

df_temp = pd.DataFrame((np.nanmean(np.array(regional_temp)[:,:,:], axis=0)))
df_temp.columns = ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW']

###############################################################################
'''
Statistic 1 (first results paragraph):
    Where most (>50%) of glacier ice exposed in Southwest and West Greenland
'''
df_count.iloc[7:12].sum() / df_count.sum()
df_count.iloc[10:14].sum() / df_count.sum()
###############################################################################


###############################################################################
'''
Statistic 2 (second results paragraph):
    Glacier ice exposure dates at 700-800 m a.s.l.
'''
df.iloc[8].sort_values()


fig, ax = plt.subplots()
ax.scatter(df_temp.iloc[8], df.iloc[8])
n = ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW']
for i, txt in enumerate(n):
    ax.annotate(txt, (df_temp.iloc[8].iloc[i], df.iloc[8].iloc[i]))
###############################################################################


























