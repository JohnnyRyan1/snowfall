#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" 

DESCRIPTION

1. Sample timing of glacier ice exposure, snow depth, and temperature from MAR.

"""

#%%

# Import modules
import xarray as xr
import netCDF4
import numpy as np
import pandas as pd

# Define user
user = 'jryan4'

# Define path
path = '/Users/' + user + '/Dropbox (University of Oregon)/research/snowfall/data/'
path_alt = '/Users/' + user + '/Dropbox (University of Oregon)/published/clouds/data/'

# Import regions
regions_file = netCDF4.Dataset(path_alt + 'temp_albedo_summer_climatologies.nc')
regions = regions_file.variables['regions'][:].filled()

# Import ISMIP 1 km grid
ismip_1km = xr.open_dataset(path + 'masks/ismip-gimp.nc')
ismip_x, ismip_y = np.meshgrid(ismip_1km['x'].values, ismip_1km['y'].values)

# Define years
years = np.arange(2001, 2022)


#%%

for year in years:
    
    print('Processing... %.0f' %year)
    
    # Import MAR grids
    mar = xr.open_dataset(path + 'mar_resample/mar_' + str(year) + '.nc')
    mar_prev = xr.open_dataset(path + 'mar_resample/mar_' + str(year - 1) + '.nc')
    
    # Import MODIS grids
    modis = xr.open_dataset(path + 'modis_exposure_v2/MOD10A1_albedo_dates_' + str(year) + '.nc')
    
    # Produce glacier ice masks
    ice_mask_first_55 = modis['first_55'].values.astype(float)
    is_mask = (ismip_1km['GIMP'].values == 0)
    ice_mask_first_55[is_mask] = np.nan
    
    ice_mask_albedo = modis['albedo'].values.astype(float)
    ice_mask_albedo[is_mask] = np.nan
    
    ice_mask_elevation = ismip_1km['SRF'].values
    
    # Find MAR grid cells that overlap with ISMIP grid
    first_55_count, first_55_median, first_55_mean = [], [], []
    albedo_median, albedo_mean = [], [] 
    elevation, region = [], []
    jun_pdd, jul_pdd, spr_temp, jun_temp = [], [], [], []
    max_snow, max_snow_day, snow_sum = [], [], []
    grid_cell_i, grid_cell_j, lon, lat = [], [], [], []
    
    # Loop over every MAR grid cell
    for i in range(mar['x'].shape[0] - 1):
        for j in range(mar['y'].shape[0] - 1):
            # Get coordinates corners
            min_y, min_x, max_y, max_x = mar['y'].values[j], mar['x'].values[i],\
                                         mar['y'].values[j+1], mar['x'].values[i+1]
    
            # Get values from MODIS
            array = ((ismip_y >= min_y) & (ismip_y < max_y) &\
                     (ismip_x >= min_x) & (ismip_x < max_x)).nonzero()[0].shape[0]
            
            if array > 0:
                
                mask = ((ismip_y >= min_y) & (ismip_y < max_y) &\
                         (ismip_x >= min_x) & (ismip_x < max_x)).nonzero()
                
                # Count valid pixels
                count = np.isfinite(ice_mask_first_55[mask]).sum()
                
                if (count > 0) & (np.nansum(ice_mask_first_55[mask]) > 0):
                    #print('stoppping...')
                    #print(i,j)
                    if np.isfinite(mar['t2m'][:, j, i].values).sum() > 0:
                        if np.isfinite(mar['sd2'][:, j, i].values).sum() > 0:
                                           
                            # June positive degree day
                            pdd = mar['t2m'][:, j, i]
                            pdd[pdd < 0] = 0
                                                        
                            jun_pdd.append(np.mean(pdd[151:181]).values)
                            
                            # July positive degree day
                            jul_pdd.append(np.mean(pdd[181:212]).values)
                            
                            # Mean Spring air temperature
                            t = mar['t2m'][:, j, i]
                            spr_temp.append(np.mean(t[60:151]).values)
                            
                            # Mean June air temperature
                            jun_temp.append(np.mean(t[151:181]).values)
                            
                            # Max snow depth, max snow depth day
                            max_snow.append(mar['sd2'][:, j, i].values[mar['sd2'][:, j, i].argmax().values])
                            max_snow_day.append(mar['time'][[mar['t2m'][:, j, i].argmax().values]].values[0])
                            
                            # Snowfall from Oct 1 to May 31
                            snow_date = mar['sd2'][:, j, i].argmax().values
                            sf_curr = mar['sf'][:,j,i][0:snow_date].sum().values
                            sf_prev = mar_prev['sf'][:, j, i][273:].sum().values
                            snow_sum.append(sf_curr + sf_prev)
                            
                            # Grid cell
                            grid_cell_i.append(i)
                            grid_cell_j.append(j)
                            lon.append(np.mean((min_x, max_x)))
                            lat.append(np.mean((min_y, max_y)))
                            
                            # Count of valid pixels
                            first_55_count.append(np.isfinite(ice_mask_first_55[mask]).sum())
                            
                            # Mean and median of valid pixels
                            first_55_median.append(np.nanmedian(ice_mask_first_55[mask][ice_mask_first_55[mask] != 0]))
                            first_55_mean.append(np.nanmean(ice_mask_first_55[mask][ice_mask_first_55[mask] != 0]))
                                                
                            # Albedo
                            albedo_median.append(np.nanmedian(ice_mask_albedo[mask][ice_mask_first_55[mask] != 0]))
                            albedo_mean.append(np.nanmean(ice_mask_albedo[mask][ice_mask_first_55[mask] != 0]))
                            
                            # Elevation and region
                            ice_mask_ele = ice_mask_elevation[mask][np.isfinite(ice_mask_first_55[mask])]
                            elevation.append(np.nanmedian(ice_mask_ele[ice_mask_ele != 0]))
                            region.append(np.nanmedian(regions[mask][np.isfinite(ice_mask_first_55[mask])]))
                        else:
                            pass
                    else:
                        pass
            
                else:
                    pass
            
            else:
                pass


    # Put in DataFrame
    df = pd.DataFrame(list(zip(lon, lat, first_55_count, first_55_median, first_55_mean,
                               albedo_median, albedo_mean, elevation, region, 
                               max_snow, max_snow_day, jun_pdd, jul_pdd, spr_temp, jun_temp, snow_sum,
                               grid_cell_i, grid_cell_j)))
            
    df.columns = ['lon', 'lat', 'first_55_count', 'first_55_median', 'first_55_mean', 
                  'albedo_median', 'albedo_mean', 'elevation', 'region', 
                  'max_snow', 'max_snow_day', 'jun_pdd', 'jul_pdd', 
                  'spr_temp', 'jun_temp', 'snow_sum',
                  'grid_cell_i', 'grid_cell_j']
    
    # Re-index
    df.set_index(['grid_cell_i', 'grid_cell_j'], inplace=True)
    
    # Save as csv
    df.to_csv(path + 'mar_results/results_' + str(year) + '.csv')






















