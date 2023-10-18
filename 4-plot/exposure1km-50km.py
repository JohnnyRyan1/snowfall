#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" 

DESCRIPTION

1. Sample timing of glacier ice exposure to show variation in each 50 x 50 km grid cell.

"""

#%%

# Import modules
import xarray as xr
import netCDF4
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

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

# Read stats DataFrame
all_stats = pd.read_csv(path + 'all_stats.csv')

#%%

for year in years:
    
    print('Processing... %.0f' %year)
    
    # Import MERRA-2 grids
    merra_curr = xr.open_dataset(path + 'merra_resample/merra_' + str(year) + '.nc')
    merra_prev = xr.open_dataset(path + 'merra_resample/merra_' + str(year - 1) + '.nc')

    # Import MODIS grids
    modis = xr.open_dataset(path + 'modis_exposure_v2/MOD10A1_albedo_dates_' + str(year) + '.nc')
    
    # Produce glacier ice masks
    ice_mask_first_55 = modis['first_55'].values.astype(float)
    is_mask = (ismip_1km['GIMP'].values == 0)
    ice_mask_first_55[is_mask] = np.nan


#%%
year = 2016

# Loop over every MERRA grid cell
for cell in range(all_stats.shape[0]):
        
    i = all_stats['grid_cell_i'].iloc[cell]
    j = all_stats['grid_cell_j'].iloc[cell]
    
    # Get coordinates corners
    min_y, min_x = merra_curr['y'].values[j], merra_curr['x'].values[i]
    max_y, max_x = merra_curr['y'].values[j+1], merra_curr['x'].values[i+1]
        
    # Get values from MODIS
    array = ((ismip_y >= min_y) & (ismip_y < max_y) &\
             (ismip_x >= min_x) & (ismip_x < max_x)).nonzero()[0].shape[0]
        
            
    mask = ((ismip_y >= min_y) & (ismip_y < max_y) &\
             (ismip_x >= min_x) & (ismip_x < max_x)).nonzero()
        
    
    data = ice_mask_first_55[mask]
    data[data == 0] = np.nan
  
    plt.hist(data, bins=30, edgecolor='k', zorder=2)
    plt.grid(zorder=1, ls='dashed')
    plt.axvline(x=np.nanmedian(data), zorder=3, color='k', lw=2, ls='dashed')
    plt.xlabel('Timing of exposure (DOY)')
    plt.ylabel('Frequency')
    plt.savefig(path + 'exposure-figures/fig_' + str(cell) + '_' + str(i) + '_' + str(j) + '.png')
    plt.close()

#%%

cells = [115, 21, 198, 101, 64, 20]

# N, N, NE, SW, CW, NW

year = 2016

data_list= []
# Loop over every MERRA grid cell
for cell in cells:
        
    i = all_stats['grid_cell_i'].iloc[cell]
    j = all_stats['grid_cell_j'].iloc[cell]
    
    # Get coordinates corners
    min_y, min_x = merra_curr['y'].values[j], merra_curr['x'].values[i]
    max_y, max_x = merra_curr['y'].values[j+1], merra_curr['x'].values[i+1]
        
    # Get values from MODIS
    array = ((ismip_y >= min_y) & (ismip_y < max_y) &\
             (ismip_x >= min_x) & (ismip_x < max_x)).nonzero()[0].shape[0]
        
            
    mask = ((ismip_y >= min_y) & (ismip_y < max_y) &\
             (ismip_x >= min_x) & (ismip_x < max_x)).nonzero()
        
    
    data = ice_mask_first_55[mask]
    data[data == 0] = np.nan
    
    data_list.append(data)

#%%
savepath = '/Users/jryan4/Dropbox (University of Oregon)/research/snowfall/re-revision/'


fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(nrows=2, ncols=3, 
                                                       figsize=(11,5), 
                                                       layout='constrained', 
                                                       sharex=True)

ax1.hist(data_list[0], bins=25, edgecolor='k', zorder=2)
ax1.axvline(x=np.nanmedian(data_list[0]), zorder=3, color='k', lw=2, ls='dashed')
ax2.hist(data_list[1], bins=15, edgecolor='k', zorder=2)
ax2.axvline(x=np.nanmedian(data_list[1]), zorder=3, color='k', lw=2, ls='dashed')
ax3.hist(data_list[2], bins=25, edgecolor='k', zorder=2)
ax3.axvline(x=np.nanmedian(data_list[2]), zorder=3, color='k', lw=2, ls='dashed')
ax4.hist(data_list[3], bins=25, edgecolor='k', zorder=2)
ax4.axvline(x=np.nanmedian(data_list[3]), zorder=3, color='k', lw=2, ls='dashed')
ax5.hist(data_list[4], bins=30, edgecolor='k', zorder=2)
ax5.axvline(x=np.nanmedian(data_list[4]), zorder=3, color='k', lw=2, ls='dashed')
ax6.hist(data_list[5], bins=30, edgecolor='k', zorder=2)
ax6.axvline(x=np.nanmedian(data_list[5]), zorder=3, color='k', lw=2, ls='dashed')

for ax in [ax1, ax2, ax3, ax4, ax5, ax6]:
    ax.grid(zorder=1, ls='dashed')

ax4.set_xlabel('Timing of exposure (DOY)', fontsize=12)
ax5.set_xlabel('Timing of exposure (DOY)', fontsize=12)
ax6.set_xlabel('Timing of exposure (DOY)', fontsize=12)
ax1.set_ylabel('Frequency', fontsize=12)
ax4.set_ylabel('Frequency', fontsize=12)

ax1.text(0.03, 0.89, "a", fontsize=20, transform=ax1.transAxes)
ax2.text(0.03, 0.87, "b", fontsize=20, transform=ax2.transAxes)
ax3.text(0.03, 0.89, "c", fontsize=20, transform=ax3.transAxes)
ax4.text(0.03, 0.87, "d", fontsize=20, transform=ax4.transAxes)
ax5.text(0.87, 0.89, "e", fontsize=20, transform=ax5.transAxes)
ax6.text(0.03, 0.87, "f", fontsize=20, transform=ax6.transAxes)

# =============================================================================
# ax1.text(0.03, 0.72, "N", fontsize=16, transform=ax1.transAxes)
# ax2.text(0.03, 0.70, "N", fontsize=16, transform=ax2.transAxes)
# ax3.text(0.03, 0.72, "NE", fontsize=16, transform=ax3.transAxes)
# ax4.text(0.03, 0.70, "SW", fontsize=16, transform=ax4.transAxes)
# ax5.text(0.82, 0.72, "CW", fontsize=16, transform=ax5.transAxes)
# ax6.text(0.03, 0.70, "NW", fontsize=16, transform=ax6.transAxes)
# 
# =============================================================================
fig.savefig(savepath + 'fig_sx_exposure-plots.png', dpi=300)

















