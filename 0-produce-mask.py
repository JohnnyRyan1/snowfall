#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" 

DESCRIPTION

1. Make masks for project 

"""

#%%

# Import modules
import xarray as xr
import geopandas as gpd
import rasterio
from rasterio import mask as msk
from rasterio.transform import Affine
import numpy as np
from skimage.morphology import binary_dilation
from skimage.measure import label, regionprops

# Define user
user = 'johnnyryan'

# Define path
path = '/Users/' + user + '/Dropbox (University of Oregon)/research/snowfall/data/'

#%%

# Read 
mask = gpd.read_file(path + 'masks/glims_polygons.shp')

# Explode
mask = mask.explode(index_parts=True)

# Reproject
mask = mask.to_crs('EPSG:3413')

#%%

# Read ISMIP
ismip = xr.open_dataset(path + 'masks/1km-ISMIP6.nc')

#%%

# Read glacier ice exposure
ice = xr.open_dataset(path + 'modis_exposure_v2/MOD10A1_albedo_dates_2010.nc')


#%%

xres = (ismip['x'].values[-1] - ismip['x'].values[0]) / len(ismip['x'].values)
yres = (ismip['y'].values[-1] - ismip['y'].values[0]) / len(ismip['y'].values)

transform = Affine.translation(ismip['x'].values[0] - xres / 2, 
                               ismip['y'].values[0] - yres / 2) * Affine.scale(xres, yres)

#%%

with rasterio.open(
        path + "masks/ismip-ice-mask.tif",
        mode="w",
        driver="GTiff",
        height=ismip['MSK'].values.shape[0],
        width=ismip['MSK'].values.shape[1],
        count=1,
        dtype=ismip['MSK'].values.dtype,
        crs="EPSG:3413",
        transform=transform,
) as new_dataset:
        new_dataset.write(ismip['MSK'].values, 1)
        
#%%

with rasterio.open(
        path + "masks/modis-ice-exposure.tif",
        mode="w",
        driver="GTiff",
        height=ismip['MSK'].values.shape[0],
        width=ismip['MSK'].values.shape[1],
        count=1,
        dtype=ismip['MSK'].values.dtype,
        crs="EPSG:3413",
        transform=transform,
) as new_dataset:
        new_dataset.write(ice['first_55'].values, 1)
        
#%%

ismip_src = rasterio.open(path + "masks/ismip-ice-mask.tif")
ismip_mask = ismip_src.read(1)

gimp_src = rasterio.open(path + 'masks/GimpIceMask_1km_2015_v1.2.tif')
gimp_mask = gimp_src.read(1)


#%%

# Clip
clipped_array, clipped_transform = msk.mask(gimp_src, mask.geometry, 
                                            crop=False, all_touched=True)

#%%

# Dilate
dilate1 = binary_dilation(np.flipud(clipped_array[0,:,:]))
dilate2 = binary_dilation(dilate1)

#%%

with rasterio.open(
        path + "masks/ismip-pgic.tif",
        mode="w",
        driver="GTiff",
        height=ismip['MSK'].values.shape[0],
        width=ismip['MSK'].values.shape[1],
        count=1,
        dtype=ismip['MSK'].values.dtype,
        crs="EPSG:3413",
        transform=transform,
) as new_dataset:
        new_dataset.write(dilate2, 1)


#%%

# Any GIMP mask pixel that is also below a dilated PGIC
pgic = (np.flipud(gimp_mask) == 1) & (dilate2 == True)
ice_sheet = (np.flipud(gimp_mask) == 1) & (dilate2 == False)


#%%
is_label = label(ice_sheet)

pixel_count = []
pixel_label = []
for l in np.unique(is_label)[0:-1]:
    pixel_count.append(regionprops(is_label)[l].num_pixels)
    pixel_label.append(regionprops(is_label)[l].label)
    
pixel_label = np.array(pixel_label)
big_labels = pixel_label[np.array(pixel_count) > 100]

ice_sheet_filter = np.isin(is_label, big_labels)

#%%

with rasterio.open(
        path + "masks/ismip-ice.tif",
        mode="w",
        driver="GTiff",
        height=ismip['MSK'].values.shape[0],
        width=ismip['MSK'].values.shape[1],
        count=1,
        dtype=ismip['MSK'].values.dtype,
        crs="EPSG:3413",
        transform=transform,
) as new_dataset:
        new_dataset.write(ice_sheet_filter, 1)
        
        #%%

ismip['PGIC'] = (['Y', 'X'],  pgic)
ismip['GIMP'] = (['Y', 'X'], ice_sheet_filter)


#%%
ismip.to_netcdf(path + 'masks/ismip-gimp.nc')















