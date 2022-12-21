""" 

DESCRIPTION

1. Sample glacier ice exposure, snowfall, and temperature on MERRA-2 grid

"""

# Import modules
import xarray as xr
import netCDF4
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Define path
path = '/Users/jryan4/Dropbox (University of Oregon)/research/snowfall/data/'

# Import regions
regions_file = netCDF4.Dataset('/Users/jryan4/Dropbox (University of Oregon)/research/clouds/data/temp_albedo_summer_climatologies.nc')
regions = regions_file.variables['regions'][:].filled()

# Import ISMIP 1 km grid
ismip_1km = xr.open_dataset(path + 'masks/1km-ISMIP6.nc')

# Define years
years = np.arange(2001, 2022)

for year in years:
    
    print('Processing... %.0f' %year)
    
    # Import MERRA-2 grids
    merra_sf = xr.open_dataset(path + 'merra_snowfall/sf_' + str(year) + '.nc')
    merra_t2m = xr.open_dataset(path + 'merra_t2m/t2m_'  + str(year) + '.nc')
    
    # Import MODIS grids
    modis = xr.open_dataset(path + 'modis_exposure_v2/MOD10A1_albedo_dates_' + str(year) + '.nc')
    
    # Produce glacier ice masks
    ice_mask_first_55 = modis['first_55'].values.astype(float)
    ice_mask_first_55[ismip_1km['ICE'] == 0] = np.nan
    
    ice_mask_first_60 = modis['first_60'].values.astype(float)
    ice_mask_first_60[ismip_1km['ICE'] == 0] = np.nan
    
    ice_mask_last_55 = modis['last_55'].values.astype(float)
    ice_mask_last_55[ismip_1km['ICE'] == 0] = np.nan
    
    ice_mask_last_60 = modis['last_60'].values.astype(float)
    ice_mask_last_60[ismip_1km['ICE'] == 0] = np.nan
    
    ice_mask_albedo = modis['albedo'].values.astype(float)
    ice_mask_albedo[ismip_1km['ICE'] == 0] = np.nan
    
    ice_mask_elevation = ismip_1km['SRF'].values
    
    # Find MERRA-2 grid cells that overlap with ISMIP grid
    first_55_count, first_60_count, last_55_count, last_60_count = [], [], [], []
    first_55_median, first_55_mean, first_60_median, first_60_mean = [], [], [], []
    last_55_median, last_55_mean, last_60_median, last_60_mean = [], [], [], []
    albedo_median, albedo_mean = [], [] 
    elevation, region = [], []
    snowfall_mean, snowfall_exposure, snowfall_coverage = [], [], []
    temp_mean, temp_exposure, temp_coverage = [], [], []
    grid_cell_i, grid_cell_j, lon, lat = [], [], [], []
    
    # Loop over every MERRA grid cell
    for i in range(merra_sf['lon'].shape[0] - 1):
        for j in range(merra_sf['lat'].shape[0] - 1):
            # Get coordinates corners
            min_y, min_x, max_y, max_x = merra_sf['latitude'].values[j], merra_sf['longitude'].values[i],\
                                         merra_sf['latitude'].values[j+1], merra_sf['longitude'].values[i+1]
    
            # Get values from MODIS
            array = ((modis.latitude >= min_y) & (modis.latitude < max_y) &\
                     (modis.longitude >= min_x) & (modis.longitude < max_x)).values.nonzero()[0].shape[0]
            
            if array > 0:
                
                mask = ((modis.latitude >= min_y) & (modis.latitude < max_y) &\
                         (modis.longitude >= min_x) & (modis.longitude < max_x)).values.nonzero()
                
                # Count valid pixels
                count = np.isfinite(ice_mask_first_55[mask]).sum()
                
                if (count > 0) & (np.nansum(ice_mask_first_55[mask]) > 0):
                
                    
                    # Count of valid pixels
                    first_55_count.append(np.isfinite(ice_mask_first_55[mask]).sum())
                    first_60_count.append(np.isfinite(ice_mask_first_60[mask]).sum())
                    last_55_count.append(np.isfinite(ice_mask_last_55[mask]).sum())
                    last_60_count.append(np.isfinite(ice_mask_last_60[mask]).sum())
                    
                    # Mean and median of valid pixels
                    first_55_median.append(np.nanmedian(ice_mask_first_55[mask][ice_mask_first_55[mask] != 0]))
                    first_55_mean.append(np.nanmean(ice_mask_first_55[mask][ice_mask_first_55[mask] != 0]))
                    first_60_median.append(np.nanmedian(ice_mask_first_60[mask][ice_mask_first_60[mask] != 0]))
                    first_60_mean.append(np.nanmean(ice_mask_first_60[mask][ice_mask_first_60[mask] != 0]))
                    
                    last_55_median.append(np.nanmedian(ice_mask_last_55[mask][ice_mask_last_55[mask] != 0]))
                    last_55_mean.append(np.nanmean(ice_mask_last_55[mask][ice_mask_last_55[mask] != 0]))
                    last_60_median.append(np.nanmedian(ice_mask_last_60[mask][ice_mask_last_60[mask] != 0]))
                    last_60_mean.append(np.nanmean(ice_mask_last_60[mask][ice_mask_last_60[mask] != 0]))
                    
                    # Albedo
                    albedo_median.append(np.nanmedian(ice_mask_albedo[mask][ice_mask_last_55[mask] != 0]))
                    albedo_mean.append(np.nanmean(ice_mask_albedo[mask][ice_mask_last_55[mask] != 0]))
                    
                    # Elevation and region
                    ice_mask_ele = ice_mask_elevation[mask][np.isfinite(ice_mask_first_55[mask])]
                    elevation.append(np.nanmedian(ice_mask_ele[ice_mask_ele != 0]))
                    region.append(np.nanmedian(regions[mask][np.isfinite(ice_mask_first_55[mask])]))
                    
                    # Get day of year index
                    idx_exp = int(np.nanmedian(ice_mask_first_55[mask][ice_mask_first_55[mask] != 0]) - 30 - 121)
                    if idx_exp < 0:
                        idx_exp = 0
                    
                    idx_cov = int(np.nanmedian(ice_mask_last_55[mask][ice_mask_last_55[mask] != 0]) - 121)
                    if idx_cov > 122:
                        idx_cov = 122
                                           
                    
                    # Mean summer snowfall, snowfall before July 1, snowfall 30 days before exposure, and 30 days before coverage
                    snowfall_mean.append(np.nanmean(merra_sf['sf'][31:123, j, i]))
                    snowfall_exposure.append(np.nanmean(merra_sf['sf'][idx_exp:idx_exp+30, j, i]))
                    snowfall_coverage.append(np.nanmean(merra_sf['sf'][idx_cov:idx_cov+30, j, i]))
                    
                    # Mean summer temperature, temperature before July 1, temperature 30 days before exposure, and 30 days before coverage
                    temp_mean.append(np.nanmean(merra_t2m['t2m'][31:123, j, i]))
                    temp_exposure.append(np.nanmean(merra_t2m['t2m'][idx_exp:idx_exp+30, j, i]))
                    temp_coverage.append(np.nanmean(merra_t2m['t2m'][idx_cov:idx_cov+30, j, i]))
                    
                    # Grid cell
                    grid_cell_i.append(i)
                    grid_cell_j.append(j)
                    lon.append(np.mean((min_x, max_x)))
                    lat.append(np.mean((min_y, max_y)))
            
                else:
                    pass
            
            else:
                pass
                
    
    # Put in DataFrame
    df = pd.DataFrame(list(zip(lon, lat, first_55_count, first_60_count, last_55_count, last_60_count,
                               first_55_median, first_60_median, first_60_mean, first_60_mean,
                               last_55_median, last_60_median, last_60_mean, last_60_mean,
                               albedo_median, albedo_mean, elevation, region, 
                               snowfall_mean, snowfall_exposure, snowfall_coverage,
                               temp_mean, temp_exposure, temp_coverage,
                               grid_cell_i, grid_cell_j)))
            
    df.columns = ['lon', 'lat', 'first_55_count', 'first_60_count', 'last_55_count', 'last_60_count',
                  'first_55_median', 'first_60_median', 'first_60_mean', 'first_60_mean',
                  'last_55_median', 'last_60_median', 'last_60_mean', 'last_60_mean',
                  'albedo_median', 'albedo_mean', 'elevation', 'region', 
                  'snowfall_mean', 'snowfall_exposure', 'snowfall_coverage',
                  'temp_mean', 'temp_exposure', 'temp_coverage',
                  'grid_cell_i', 'grid_cell_j']
    
    # Re-index
    df.set_index(['grid_cell_i', 'grid_cell_j'], inplace=True)
    
    # Save as csv
    df.to_csv(path + 'results/results_' + str(year) + '.csv')





















