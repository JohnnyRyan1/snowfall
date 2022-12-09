import xarray as xr
import numpy as np
from pyresample.bucket import BucketResampler
from pyresample import create_area_def

# Define base path
path = '/Users/jryan4/Dropbox (University of Oregon)/research/snowfall/'

# Define MERRA data
merra_file = xr.open_dataset(path + 'data/merra_sample/MERRA2_200.tavg1_2d_int_Nx.20000101.SUB.nc')

# Read ISMIP data
ismip = xr.open_dataset(path + 'data/masks/1km-ISMIP6.nc')

# Meshgrid lat/lons
merra_mesh_lon, merra_mesh_lat = np.meshgrid(merra_file['lon'], merra_file['lat'])

# Create area definition
dy = 0.5
dx = 0.625
target_def = create_area_def('merra2',
                           {'proj': 'longlat', 'datum': 'WGS84'},
                           area_extent=[merra_file['lon'].data.min()-dx/2, merra_file['lat'].data.min()-dy/2, 
                                        merra_file['lon'].data.max()+dx/2, merra_file['lat'].data.max()+dy/2],
                           resolution=(dx, dy),
                           units='degrees',
                           description='Greenland 0.625 x 0.5 degree lat-lon grid')

# Resample
resampler = BucketResampler(target_def, ismip['lon'].data, ismip['lat'].data)



















