#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" 

DESCRIPTION

1. Plot showing relationship between snow depth, cumulative temperature, 
and timing of glacier ice exposure.

"""

# Import modules
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

# Define path
path = '/Users/jryan4/Dropbox (University of Oregon)/research/snowfall/data/'

# Define save path
savepath = '/Users/jryan4/Dropbox (University of Oregon)/research/snowfall/manuscript/figures/'

# Import datasets
modis = xr.open_dataset(path + 'modis_exposure_v2/MOD10A1_albedo_dates_2010.nc', engine='netcdf4')
merra_sd = xr.open_dataset(path + 'merra_snowdepth/sd_2010.nc')
merra_t2m = xr.open_dataset(path + 'merra_t2m/t2m_2010.nc') 

# Import ISMIP 1 km grid
ismip_1km = xr.open_dataset(path + 'masks/1km-ISMIP6.nc')

# Produce glacier ice masks
ice_mask_first_55 = modis['first_55'].values.astype(float)
ice_mask_first_55[ismip_1km['ICE'] == 0] = np.nan
ice_mask_albedo = modis['albedo'].values.astype(float)
ice_mask_albedo[ismip_1km['ICE'] == 0] = np.nan
ice_mask_elevation = ismip_1km['SRF'].values

# Find median indices
medians = np.argwhere(ice_mask_first_55 == 183)

# Get single pixel time-series
pixel = 231
#print(medians[pixel][0], medians[pixel][1])
expo_test = ice_mask_first_55[medians[pixel][0], medians[pixel][1]]
#plt.imshow(ice_mask_first_55)
#plt.scatter(medians[pixel,1], medians[pixel,0])

# Find which MERRA grid cell contains this grid cell
modis_point = [modis['latitude'].values[medians[pixel][0], medians[pixel][1]],
               modis['longitude'].values[medians[pixel][0], medians[pixel][1]]]

abslat = np.abs(merra_sd.latitude - modis_point[0])
abslon = np.abs(merra_sd.longitude - modis_point[1])
c = np.maximum(abslon, abslat)
([xloc], [yloc]) = np.where(c == np.min(c))

# Convert temperature to positive degree days
pdd = merra_t2m['t2m'][:, yloc, xloc] - 273
pdd[pdd < 0] = 0

# Compute timing of max snow depth
max_sd = merra_sd['t2m'][:, yloc, xloc].values[merra_sd['t2m'][:, yloc, xloc].argmax().values]
max_sd_day = merra_sd['time'].dt.dayofyear[[merra_sd['t2m'][:, yloc, xloc].argmax().values]].values[0]

###############################################################################

# Plot figure
fig, (ax1) = plt.subplots(nrows=1, ncols=1, figsize=(8, 4), layout='constrained')

ax1.plot(merra_sd['time'].dt.dayofyear, merra_sd['t2m'][:, yloc, xloc], color='#0000a7', 
         lw=2, alpha=0.8, zorder=3)
ax1.fill_between(merra_sd['time'].dt.dayofyear, merra_sd['t2m'][:, yloc, xloc], 
                 where = (merra_sd['time'].dt.dayofyear > 58) & (merra_sd['time'].dt.dayofyear <= max_sd_day), 
                 alpha=0.4, color='#0000a7', zorder=2)
ax1.axvline(x=expo_test, ls='dashed', lw=2, color='k')

ax2 = ax1.twinx()
ax2.plot(merra_t2m['time'].dt.dayofyear, np.cumsum(pdd), lw=2, color='#c1272d', 
         alpha=0.8, zorder=3)
ax2.fill_between(merra_t2m['time'].dt.dayofyear, np.cumsum(pdd), 
                 where= (merra_t2m['time'].dt.dayofyear > 151) & (merra_t2m['time'].dt.dayofyear < 183), 
                 alpha=0.4, color='#c1272d', zorder=2)



ax1.set_ylabel('Snow depth (m)', fontsize=14)
ax2.set_ylabel('Cumulative PDD (K)', fontsize=14)
ax1.set_xlabel('Time (Day of year)', fontsize=14)

ax1.set_ylim(0, 1.3)
ax2.set_ylim(0, 250)
ax1.set_xlim(60, 274)
ax1.grid(ls='dashed', lw=1, zorder=1)
ax1.tick_params(axis='both', which='major', labelsize=14)
ax2.tick_params(axis='both', which='major', labelsize=14)

#fig.savefig(savepath + 'fig_sx_snow_depth_vs_temp.pdf')

###############################################################################

























