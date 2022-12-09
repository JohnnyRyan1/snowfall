""" 

DESCRIPTION

1. Sample glacier ice exposure, snowfall, and temperature on MERRA-2 grid

"""

# Import modules
import xarray as xr
import netCDF4
import numpy as np
import glob
from pyresample import geometry, utils, image

# Define path
path = '/Users/jryan4/Dropbox (University of Oregon)/research/snowfall/data/'

# Define destination to save
dest_1km = path + 'modis_exposure_v2/'

# Import MERRA-2 grid
merra = xr.open_dataset(path + 'merra_snowfall_monthly/MERRA2_400.tavgM_2d_int_Nx.201706.SUB.nc')

# Import ISMIP 1 km grid
ismip_1km = netCDF4.Dataset(path + 'masks/1km-ISMIP6.nc','r')

# Import MODIS grids



