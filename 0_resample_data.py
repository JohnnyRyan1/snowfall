#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
DESCRIPTION

1. Reproject MERRA-2 data

"""

#%%

# Import modules
import xarray as xr
from pyresample import kd_tree, geometry
import numpy as np
from scipy import ndimage
import pandas as pd
from datetime import datetime
import os

# Define user
user = 'johnnyryan'

# Define path
path = '/Users/' + user + '/Dropbox (University of Oregon)/research/snowfall/data/'
dest = path + 'merra_resample/'

#%%

# Import ISMIP 1 km grid
ismip_1km = xr.open_dataset(path + 'masks/ismip-gimp.nc')

# Covert ISMIP to 50 km grid cells
lons_50km = ndimage.zoom(ismip_1km['lon'].values, 0.02)
lats_50km = ndimage.zoom(ismip_1km['lat'].values, 0.02)
x_50km = ndimage.zoom(ismip_1km['x'].values, 0.02)
y_50km = ndimage.zoom(ismip_1km['y'].values, 0.02)

# Define target projection
target_def = geometry.SwathDefinition(lons=lons_50km, lats=lats_50km)

# Define years
years = np.arange(1980, 2022)

#%%

"""
Resample MERRA-2 1980 to 2021

"""

for year in years:
    
    # Read files
    merra_t2m = xr.open_dataset(path + 'merra_t2m/t2m_'  + str(year) + '.nc')
    merra_sf = xr.open_dataset(path + 'merra_snowfall/sf_'  + str(year) + '.nc')
    
    # Extract lat/lon
    merra_lon, merra_lat = np.meshgrid(merra_t2m['longitude'].values, merra_t2m['latitude'].values)
    
    # Define source projection
    source_def = geometry.SwathDefinition(lons=merra_lon, lats=merra_lat)
        
    # Reproject
    t2m = kd_tree.resample_nearest(source_def, np.rollaxis(merra_t2m['t2m'].values, 0, 3), 
                                      target_def, radius_of_influence=50000)
    
    # Convert zeros to NaNs
    t2m[t2m == 0] = np.nan
    
    # Roll axis back
    t2m = np.rollaxis(t2m, 2, 0)
    
    # Reproject
    sf = kd_tree.resample_nearest(source_def, np.rollaxis(merra_sf['sf'].values, 0, 3), 
                                      target_def, radius_of_influence=50000)
    
    # Convert zeros to NaNs
    sf[sf == 0] = np.nan
    
    # Roll axis back
    sf = np.rollaxis(sf, 2, 0)

    # Save as NetCDF
    ds_data = xr.Dataset(
    data_vars={
        "t2m": (("time", "y", "x"), t2m.astype('float32')),
        "sf": (("time", "y", "x"), sf.astype('float32')),
    },

    coords={
        "time": pd.DatetimeIndex(merra_t2m['time'].values, freq='D'),
        "y": (('y',), y_50km),
        "x": (('x',), x_50km),    
    },

    attrs={
        "Produced": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "Author":'Johnny Ryan', 
        "Email":'jryan4@uoregon.edu'
    },
    )
    
    # Save
    ds_data.to_netcdf(dest + 'merra_' + str(year) + '.nc')
    
#%% 

"""
Resample MAR

"""

# Define years
years = np.arange(2000, 2022)
dest = path + 'mar_resample/'

for year in years:

    if os.path.exists(dest + 'mar_' + str(year) + '.nc'):
        print(f'Skipping...{str(year)}')
    else:
        print(f'Processing...{str(year)}')
    
        # Read files
        mar = xr.open_dataset(path + 'mar-v3-12-1/MARv3.12.1-10km-daily-ERA5-' + str(year) + '.nc')
        
        # Coarsen
        mar_50km = mar.coarsen(x=5).mean().coarsen(y=5).mean()
        
        # Get snowfall
        mar_sf = mar_50km['SF'][:,:,:].values
        mar_t2m = mar_50km['TTZ'][:,0,:,:].values
    
        # Define source projection
        source_def = geometry.SwathDefinition(lons=mar_50km['LON'].values, lats=mar_50km['LAT'].values)
    
        # Reproject 
        mar_sf_resample = kd_tree.resample_nearest(source_def, np.rollaxis(mar_sf, 0, 3), 
                                          target_def, radius_of_influence=50000)
        mar_t2m_resample = kd_tree.resample_nearest(source_def, np.rollaxis(mar_t2m, 0, 3), 
                                          target_def, radius_of_influence=50000)
        
        # Convert zeros to NaNs
        mar_sf_resample[mar_sf_resample == 0] = np.nan
        mar_t2m_resample[mar_t2m_resample == 0] = np.nan
        
        # Roll axis back
        mar_sf_resample_final = np.rollaxis(mar_sf_resample, 2, 0)
        mar_t2m_resample_final = np.rollaxis(mar_t2m_resample, 2, 0)
        
        # Save as NetCDF
        ds_data = xr.Dataset(
        data_vars={
            "t2m": (("time", "y", "x"), mar_t2m_resample_final.astype('float32')),
            "sf": (("time", "y", "x"), mar_sf_resample_final.astype('float32')),
        },
    
        coords={
            "time": pd.DatetimeIndex(mar['TIME'].values, freq='D'),
            "y": (('y',), y_50km),
            "x": (('x',), x_50km),    
        },
    
        attrs={
            "Produced": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "Author":'Johnny Ryan', 
            "Email":'jryan4@uoregon.edu'
        },
        )
        
        # Save
        ds_data.to_netcdf(dest + 'mar_' + str(year) + '.nc')






















