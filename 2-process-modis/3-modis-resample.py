#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" 

DESCRIPTION

1. Resample MODIS data to 1 km ISMIP grid.

"""

# Import modules
import netCDF4
import numpy as np
import glob
from pyresample import geometry, utils, image

# Define path
path = '/Users/jryan4/Dropbox (University of Oregon)/research/snowfall/data/'

# Define destination to save
dest_1km = path + 'modis_exposure_v2/'

# Import ISMIP 1 km grid
ismip_1km = netCDF4.Dataset(path + 'masks/1km-ISMIP6.nc','r')

# Define MODIS files
modis_files = sorted(glob.glob(path + 'modis_exposure_intermediate/*.nc'))
print('Number of files in folder = %.0f' %len(modis_files))

# Define years
years = np.arange(2001, 2022, 1)

for i in years:
    modis_list = []   
    # Get MODIS tiles
    for f in modis_files:
        if  f[-7:-3]  == str(i):
            modis_list.append(f)
    
    # Define new master grid
    master_grid_first_60 = np.zeros((7200, 7200), dtype='int16')
    master_grid_first_55 = np.zeros((7200, 7200), dtype='int16')
    master_grid_last_60 = np.zeros((7200, 7200), dtype='int16')
    master_grid_last_55 = np.zeros((7200, 7200), dtype='int16')
    master_grid_albedo = np.zeros((7200,7200), dtype='float')
    master_grid_lat = np.zeros((7200, 7200), dtype='float')
    master_grid_lon = np.zeros((7200, 7200), dtype='float')

    # Add tile to master grid
    for j in modis_list:
        if j[-13:-8] == '15v00':
            modis = netCDF4.Dataset(j, 'r')
            master_grid_first_55[0:2400,0:2400] = modis.variables['first_55'][:]
            master_grid_first_60[0:2400,0:2400] = modis.variables['first_60'][:]
            master_grid_last_55[0:2400,0:2400] = modis.variables['last_55'][:]
            master_grid_last_60[0:2400,0:2400] = modis.variables['last_60'][:]
            master_grid_albedo[0:2400,0:2400] = modis.variables['albedo'][:]
            master_grid_lat[0:2400,0:2400] = modis.variables['latitude'][:]
            master_grid_lon[0:2400,0:2400] = modis.variables['longitude'][:]
        if j[-13:-8] == '16v00':
            modis = netCDF4.Dataset(j, 'r')
            master_grid_first_55[0:2400,2400:4800] = modis.variables['first_55'][:]
            master_grid_first_60[0:2400,2400:4800] = modis.variables['first_60'][:]
            master_grid_last_55[0:2400,2400:4800] = modis.variables['last_55'][:]
            master_grid_last_60[0:2400,2400:4800] = modis.variables['last_60'][:]
            master_grid_albedo[0:2400,2400:4800] = modis.variables['albedo'][:]
            master_grid_lat[0:2400,2400:4800] = modis.variables['latitude'][:]
            master_grid_lon[0:2400,2400:4800] = modis.variables['longitude'][:]
        if j[-13:-8] == '17v00':
            modis = netCDF4.Dataset(j, 'r')
            master_grid_first_55[0:2400,4800:7200] = modis.variables['first_55'][:]
            master_grid_first_60[0:2400,4800:7200] = modis.variables['first_60'][:]
            master_grid_last_55[0:2400,4800:7200] = modis.variables['last_55'][:]
            master_grid_last_60[0:2400,4800:7200] = modis.variables['last_60'][:]
            master_grid_albedo[0:2400,4800:7200] = modis.variables['albedo'][:]
            master_grid_lat[0:2400,4800:7200] = modis.variables['latitude'][:]
            master_grid_lon[0:2400,4800:7200] = modis.variables['longitude'][:]
        if j[-13:-8] == '15v01':
            modis = netCDF4.Dataset(j, 'r')
            master_grid_first_55[2400:4800,0:2400] = modis.variables['first_55'][:]
            master_grid_first_60[2400:4800,0:2400] = modis.variables['first_60'][:]
            master_grid_last_55[2400:4800,0:2400] = modis.variables['last_55'][:]
            master_grid_last_60[2400:4800,0:2400] = modis.variables['last_60'][:]
            master_grid_albedo[2400:4800,0:2400] = modis.variables['albedo'][:]
            master_grid_lat[2400:4800,0:2400] = modis.variables['latitude'][:]
            master_grid_lon[2400:4800,0:2400] = modis.variables['longitude'][:]
        if j[-13:-8] == '16v01':
            modis = netCDF4.Dataset(j, 'r')
            master_grid_first_55[2400:4800,2400:4800] = modis.variables['first_55'][:]
            master_grid_first_60[2400:4800,2400:4800] = modis.variables['first_60'][:]
            master_grid_last_55[2400:4800,2400:4800] = modis.variables['last_55'][:]
            master_grid_last_60[2400:4800,2400:4800] = modis.variables['last_60'][:]
            master_grid_albedo[2400:4800,2400:4800] = modis.variables['albedo'][:]
            master_grid_lat[2400:4800,2400:4800] = modis.variables['latitude'][:]
            master_grid_lon[2400:4800,2400:4800] = modis.variables['longitude'][:]
        if j[-13:-8] == '17v01':
            modis = netCDF4.Dataset(j, 'r')
            master_grid_first_55[2400:4800,4800:7200] = modis.variables['first_55'][:]
            master_grid_first_60[2400:4800,4800:7200] = modis.variables['first_60'][:]
            master_grid_last_55[2400:4800,4800:7200] = modis.variables['last_55'][:]
            master_grid_last_60[2400:4800,4800:7200] = modis.variables['last_60'][:]
            master_grid_albedo[2400:4800,4800:7200] = modis.variables['albedo'][:]
            master_grid_lat[2400:4800,4800:7200] = modis.variables['latitude'][:]
            master_grid_lon[2400:4800,4800:7200] = modis.variables['longitude'][:]
        if j[-13:-8] == '15v02':
            modis = netCDF4.Dataset(j, 'r')
            master_grid_first_55[4800:7200,0:2400] = modis.variables['first_55'][:]
            master_grid_first_60[4800:7200,0:2400] = modis.variables['first_60'][:]
            master_grid_last_55[4800:7200,0:2400] = modis.variables['last_55'][:]
            master_grid_last_60[4800:7200,0:2400] = modis.variables['last_60'][:]
            master_grid_albedo[4800:7200,0:2400] = modis.variables['albedo'][:]
            master_grid_lat[4800:7200,0:2400] = modis.variables['latitude'][:]
            master_grid_lon[4800:7200,0:2400] = modis.variables['longitude'][:]
        if j[-13:-8] == '16v02':
            modis = netCDF4.Dataset(j, 'r')
            master_grid_first_55[4800:7200,2400:4800] = modis.variables['first_55'][:]
            master_grid_first_60[4800:7200,2400:4800] = modis.variables['first_60'][:]
            master_grid_last_55[4800:7200,2400:4800] = modis.variables['last_55'][:]
            master_grid_last_60[4800:7200,2400:4800] = modis.variables['last_60'][:]
            master_grid_albedo[4800:7200,2400:4800] = modis.variables['albedo'][:]
            master_grid_lat[4800:7200,2400:4800] = modis.variables['latitude'][:]
            master_grid_lon[4800:7200,2400:4800] = modis.variables['longitude'][:]
        if j[-13:-8] == '17v02':
            modis = netCDF4.Dataset(j, 'r')
            master_grid_first_55[4800:7200,4800:7200] = modis.variables['first_55'][:]
            master_grid_first_60[4800:7200,4800:7200] = modis.variables['first_60'][:]
            master_grid_last_55[4800:7200,4800:7200] = modis.variables['last_55'][:]
            master_grid_last_60[4800:7200,4800:7200] = modis.variables['last_60'][:]
            master_grid_albedo[4800:7200,4800:7200] = modis.variables['albedo'][:]
            master_grid_lat[4800:7200,4800:7200] = modis.variables['latitude'][:]
            master_grid_lon[4800:7200,4800:7200] = modis.variables['longitude'][:]

    # Get ISMIP6 lat lons
    lon_1km = ismip_1km.variables['lon'][:]
    lat_1km = ismip_1km.variables['lat'][:]
    
    # Convert 0s to NaNs so they do not interfere with resampling
    master_grid_first_55 = master_grid_first_55.astype('float')
    master_grid_first_60 = master_grid_first_60.astype('float')
    master_grid_last_55 = master_grid_last_55.astype('float')
    master_grid_last_60 = master_grid_last_60.astype('float')
    master_grid_albedo = master_grid_albedo.astype('float')
    
    master_grid_first_55[master_grid_first_55 == 0] = np.nan
    master_grid_first_60[master_grid_first_60 == 0] = np.nan
    master_grid_last_55[master_grid_last_55 == 0] = np.nan
    master_grid_last_60[master_grid_last_60 == 0] = np.nan
    master_grid_albedo[master_grid_albedo == 0] = np.nan
    
    # Define regridding conversion
    swath_def = geometry.SwathDefinition(lons=lon_1km, lats=lat_1km)
    swath_con = geometry.SwathDefinition(lons=master_grid_lon, lats=master_grid_lat)
    first_55_con = image.ImageContainer(master_grid_first_55, swath_con)
    first_60_con = image.ImageContainer(master_grid_first_60, swath_con)
    last_55_con = image.ImageContainer(master_grid_last_55, swath_con)
    last_60_con = image.ImageContainer(master_grid_last_60, swath_con)
    albedo_con = image.ImageContainer(master_grid_albedo, swath_con)
    row_indices, col_indices = utils.generate_nearest_neighbour_linesample_arrays(swath_con, swath_def, 1000)
    
    # Perform regridding
    first_55_result = first_55_con.get_array_from_linesample(row_indices, col_indices)
    first_60_result = first_60_con.get_array_from_linesample(row_indices, col_indices)
    last_55_result = last_55_con.get_array_from_linesample(row_indices, col_indices)
    last_60_result = last_60_con.get_array_from_linesample(row_indices, col_indices)
    albedo_result = albedo_con.get_array_from_linesample(row_indices, col_indices)
    
    ###############################################################################
    # Save 1 km dataset to NetCDF
    ###############################################################################
    dataset = netCDF4.Dataset(dest_1km + 'MOD10A1_albedo_dates_' + str(i) + '.nc', 
                              'w', format='NETCDF4_CLASSIC')
    print('Creating %s' %dest_1km + 'MOD10A1_albedo_dates_' + str(i) + '.nc')
    dataset.Title = "Mean albedo and dates of first and last bare ice exposure for %s from MOD10A1 product" %(str(i))
    import time
    dataset.History = "Created " + time.ctime(time.time())
    dataset.Projection = "WGS 84"
    dataset.Reference = "Ryan, J. C. et al. (unpublished)"
    dataset.Contact = "jryan4@uoregon.edu"
        
    # Create new dimensions
    lat_dim = dataset.createDimension('y', first_55_result.shape[0])
    lon_dim = dataset.createDimension('x', first_55_result.shape[1])
        
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
    Y[:] = lat_1km
    X[:] = lon_1km
    x[:] = lon_1km[0,:]
    y[:] = lat_1km[:,0]
    first_ice_55[:] = first_55_result
    first_ice_60[:] = first_60_result
    last_ice_55[:] = last_55_result
    last_ice_60[:] = last_60_result
    albedo_nc[:] = albedo_result
    
    print('Writing data to %s' %dest_1km + 'MOD10A1_albedo_dates_' + str(i) + '.nc')
        
    # Close dataset
    dataset.close()
