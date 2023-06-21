#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

1. Compute the first and last day of bare ice exposure from M0D10A1.

2. Compute mean summer albedo from MOD10A1.

"""

# Import modules
import numpy as np
import netCDF4
import xarray as xr
import os
import glob

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

# Define path
path = '/Users/jryan4/Dropbox (University of Oregon)/research/snowfall/data/'

# Define file list
modis_files = sorted(glob.glob(path + 'modis_albedo_intermediate/*.nc'))

# Define destination to save files
dest = path + 'modis_exposure_intermediate/'

def median_filter(albedo):
    
    # Get number of valid albedo observations
    valid_obs = np.count_nonzero(~np.isnan(albedo),axis=2)
    valid_obs_3d = np.repeat(valid_obs[:, :, np.newaxis], albedo.shape[2], axis=2)
    
    # Perform filter (i.e. remove if two standard deviations from mean)
    rolling_mean = albedo.rolling(z=11, min_periods=1, center=True).mean()
    rolling_std = albedo.rolling(z=11, min_periods=1, center=True).std()
    rolling_std = rolling_std * 2
    
    # Calculate difference between pixel value and rolling mean
    difference = np.abs(albedo - rolling_mean)
    
    # Mask values that are more than two standard deviations from the mean
    mask = (difference < rolling_std)
    
    # Calculate 11-day rolling median to be used as the timeseries
    rolling_median = albedo.where(mask == True).rolling(z=11, min_periods=3, center=True).median()
    
    # Replace value with less than 10% valid observations with 100
    rolling_valid = rolling_median.where((valid_obs_3d > 10), 100)

    # Calculate first and last day of bare ice exposure
    bare_ice_60 = rolling_median.where(rolling_valid < 60)
    bare_ice_55 = rolling_median.where(rolling_valid < 55)
    
    # Calculate number of bare ice days
    valid_bare_ice = np.repeat(np.count_nonzero(~np.isnan(bare_ice_55),axis=2)[:,:, np.newaxis], albedo.shape[2], axis=2)
    
    # Replace values where bare ice duration is less than 7 days with NaN
    bare_ice_60 = bare_ice_60.where((valid_bare_ice > 7), np.nan)
    bare_ice_55 = bare_ice_55.where((valid_bare_ice > 7), np.nan)

    # Replace all values with day of year and get first and last bare ice date
    first_bare_ice_60 = bare_ice_60.where(np.isnan(bare_ice_60), bare_ice_60['z']).min(axis=2).astype('int16')
    last_bare_ice_60 = bare_ice_60.where(np.isnan(bare_ice_60), bare_ice_60['z']).max(axis=2).astype('int16')

    # Replace all values with day of year and get first and last bare ice date
    first_bare_ice_55 = bare_ice_55.where(np.isnan(bare_ice_55), bare_ice_55['z']).min(axis=2).astype('int16')
    last_bare_ice_55 = bare_ice_55.where(np.isnan(bare_ice_55), bare_ice_55['z']).max(axis=2).astype('int16')

    return first_bare_ice_60, last_bare_ice_60, first_bare_ice_55, last_bare_ice_55

def save2netcdf(dest, lats, lons, first_55, first60, last_55, last_60, albedo, tile, year):
    dataset = netCDF4.Dataset(dest + 'MOD10A1_albedo_dates_' + tile + '_' + year + '.nc', 
                              'w', format='NETCDF4_CLASSIC')
    #print('Creating %s' %dest + 'MOD10A1_albedo_dates_' + tile + '_' + year + '.nc')
    dataset.Title = "Mean albedo and dates of first and last ice exposure for tile %s for %s from MOD10A1 product" %(tile, year)
    import time
    dataset.History = "Created " + time.ctime(time.time())
    dataset.Projection = "WGS 84"
    dataset.Reference = "Ryan, J. C. et al. (unpublished)"
    dataset.Contact = "jryan4@uoregon.edu"
        
    # Create new dimensions
    lat_dim = dataset.createDimension('y', first_55.shape[0])
    lon_dim = dataset.createDimension('x', first_55.shape[1])
        
    # Define variable types
    Y = dataset.createVariable('latitude', np.float64, ('y','x'))
    X = dataset.createVariable('longitude', np.float64, ('y','x'))
    
    y = dataset.createVariable('y', np.float64, ('y'))
    x = dataset.createVariable('x', np.float64, ('x'))
        
    # Define units
    Y.units = "degrees"
    X.units = "degrees"
       
    # Create the actual 3D variable
    first_ice_55 = dataset.createVariable('first_55', np.int16, ('y','x'))
    first_ice_60 = dataset.createVariable('first_60', np.int16, ('y','x'))
    last_ice_55 = dataset.createVariable('last_55', np.int16, ('y','x'))
    last_ice_60 = dataset.createVariable('last_60', np.int16, ('y','x'))
    albedo_nc = dataset.createVariable('albedo', np.float64, ('y','x'))
    
    # Write data to layers
    Y[:] = lats
    X[:] = lons
    x[:] = lons[0,:]
    y[:] = lats[:,0]
    first_ice_55[:] = first_55
    first_ice_60[:] = first60
    last_ice_55[:] = last_55
    last_ice_60[:] = last_60
    albedo_nc[:] = albedo

    #print('Writing data to %s' %dest + 'MOD10A1_albedo_dates_' + tile + '_' + year + '.nc')
        
    # Close dataset
    dataset.close()

for f in modis_files:

    # Import MODIS data
    mod = xr.open_dataset(f)
    
    if mod['day_of_year'][0].values == 152:
        
        # Get path and filename seperately 
        infilepath, infilename = os.path.split(f) 
        # Get file name without extension            
        infileshortname, extension = os.path.splitext(infilename)

        print('Processing... %s' %(f))
        if os.path.exists(dest + 'MOD10A1_albedo_dates_' + infileshortname[15:20] + '_' + infileshortname[21:25] + '.nc'):
            print('Skipping... %s' %(dest + 'MOD10A1_albedo_dates_' + infileshortname[15:20] + '_' + infileshortname[21:25] + '.nc'))
        else:    
            # Split into 100 x 100
            chunks = np.arange(0, 2500, 100)

            # Create empty arrays for new data
            first_bare_ice_60 = np.zeros(mod['albedo'].shape[0:2]).astype('int16')
            last_bare_ice_60 = np.zeros(mod['albedo'].shape[0:2]).astype('int16')
            first_bare_ice_55 = np.zeros(mod['albedo'].shape[0:2]).astype('int16')
            last_bare_ice_55 = np.zeros(mod['albedo'].shape[0:2]).astype('int16')
            albedo = np.zeros(mod['albedo'].shape[0:2]).astype('float16')

            for i in range(len(chunks) - 1):
                for j in range(len(chunks) - 1):
                    #print('%.0f and %0.f and %.0f and %0.f' %(chunks[i],chunks[i+1],chunks[j],chunks[j+1]))
                    mask = mod['mask'][chunks[i]:chunks[i+1],chunks[j]:chunks[j+1]]
                    if mask.sum() > 0:
                        #print('Mask sum = %0.f' %mask.sum())                    
                        array = mod['albedo'][chunks[i]:chunks[i+1],chunks[j]:chunks[j+1]]

                        # Convert to float so you can use NaNs
                        array = array.astype('float32')

                        # Filter out flagged values
                        array = array.where(array != 125, 25)
                        array = array.where(array < 100)
                        array = array.where(array > 20)

                        # Perform filter
                        first_60, last_60, first_55, last_55 = median_filter(array)
                        first_bare_ice_60[chunks[i]:chunks[i+1],chunks[j]:chunks[j+1]] = first_60
                        last_bare_ice_60[chunks[i]:chunks[i+1],chunks[j]:chunks[j+1]] = last_60
                        first_bare_ice_55[chunks[i]:chunks[i+1],chunks[j]:chunks[j+1]] = first_55
                        last_bare_ice_55[chunks[i]:chunks[i+1],chunks[j]:chunks[j+1]] = last_55
                        albedo[chunks[i]:chunks[i+1],chunks[j]:chunks[j+1]] = np.nanmean(array[:,:,0:92], axis=2).astype('float16')

            # Filter out if first and last bare ice exposure are the same
            first_bare_ice_60[first_bare_ice_60 == last_bare_ice_60] = 0
            first_bare_ice_55[first_bare_ice_55 == last_bare_ice_55] = 0
            last_bare_ice_60[first_bare_ice_60 == last_bare_ice_60] = 0
            last_bare_ice_55[first_bare_ice_55 == last_bare_ice_55] = 0

            # Filter out if first day of bare ice is in September or later
            first_bare_ice_60[first_bare_ice_55 > 244] = 0
            first_bare_ice_55[first_bare_ice_55 > 244] = 0
            last_bare_ice_60[first_bare_ice_55 > 244] = 0
            last_bare_ice_60[first_bare_ice_55 > 244] = 0

            first_bare_ice_60[first_bare_ice_55 < 90] = 0
            first_bare_ice_55[first_bare_ice_55 < 90] = 0
            last_bare_ice_60[first_bare_ice_55 < 90] = 0
            last_bare_ice_60[first_bare_ice_55 < 90] = 0

            # Save as NetCDF
            save2netcdf(dest, mod['latitude'], mod['longitude'], first_bare_ice_55, first_bare_ice_60,
                        last_bare_ice_55, last_bare_ice_60, albedo, infileshortname[15:20], 
                        infileshortname[21:25])

            # =========================================================================
            # # Optional export csv to test
            # first_bare_ice = first_bare_ice.astype('float')
            # first_bare_ice[first_bare_ice == 0] = np.nan
            # mask = np.ravel(~np.isnan(first_bare_ice))
            # df = pd.DataFrame(list(zip(np.ravel(lon)[mask], np.ravel(lat)[mask],
            #                                     np.ravel(first_bare_ice)[mask])))
            # df.columns = ['X', 'Y', 'Z']
            # df.to_csv('/media/johnny/RyanLab/jryan11/PolarSnow/Data/test.csv')
            # =========================================================================
    else:
        print('Data does not start on May 1')
        break
        
for f in modis_files:

    # Import MODIS data
    mod = xr.open_dataset(f)
    
    if mod['day_of_year'][0].values != 152:
        print('Bad start data... %s' %(f))
    
    else:
        pass

# Testing

# Import MODIS data
mod = xr.open_dataset(modis_files[0])

# Split into 100 x 100
chunks = np.arange(0, 2500, 100)

# Create empty arrays for new data
first_bare_ice_60 = np.zeros(mod['albedo'].shape[0:2]).astype('int16')
last_bare_ice_60 = np.zeros(mod['albedo'].shape[0:2]).astype('int16')
first_bare_ice_55 = np.zeros(mod['albedo'].shape[0:2]).astype('int16')
last_bare_ice_55 = np.zeros(mod['albedo'].shape[0:2]).astype('int16')
albedo = np.zeros(mod['albedo'].shape[0:2]).astype('float16')

i=23
j=12