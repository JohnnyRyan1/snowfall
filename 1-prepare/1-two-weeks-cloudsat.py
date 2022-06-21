"""

Collate all CloudSat tracks for a single cycle.

"""


# Import packages
import geopandas as gpd
import glob
import numpy as np
import pandas as pd
from pyhdf.SD import SD, SDC
from pyhdf import HDF, VS, V
import glob
import os

# Define filepaths
fp = '/Users/jryan4/Dropbox (University of Oregon)/research/snowfall/data/'

# Define file list
files = sorted(glob.glob(fp + '2c-snow-profile-r05/*.hdf'))

# Define ROI
roi = gpd.read_file(fp + 'masks/GRL_adm0.shp')

lon_data = []
lat_data = []

for i in files[0:400]:
    # print(i)
    # Extract lats and lons
    h = HDF.HDF(i)
    vs = h.vstart()

    xid = vs.find('Latitude')
    latid = vs.attach(xid)
    latid.setfields('Latitude')
    nrecs, _, _, _, _ = latid.inquire()
    latitude = latid.read(nRec=nrecs)
    latid.detach()
    latitude = np.asarray(latitude)

    lonid = vs.attach(vs.find('Longitude'))
    lonid.setfields('Longitude')
    nrecs, _, _, _, _ = lonid.inquire()
    longitude = lonid.read(nRec=nrecs)
    lonid.detach()
    longitude = np.asarray(longitude)

    sfid = vs.attach(vs.find('snowfall_rate_sfc'))
    sfid.setfields('snowfall_rate_sfc')
    nrecs, _, _, _, _ = sfid.inquire()
    snowfall_rate = sfid.read(nRec=nrecs)
    sfid.detach()
    snowfall_rate = np.asarray(snowfall_rate)
    snowfall_rate[snowfall_rate == -999] = np.nan

    tid = vs.attach(vs.find('TAI_start'))
    tid.setfields('TAI_start')
    nrecs, _, _, _, _ = tid.inquire()
    time = tid.read(nRec=nrecs)
    tid.detach()
    # seconds since 00:00:00 Jan 1 1993
    time = np.asarray(time)

    pid = vs.attach(vs.find('Profile_time'))
    pid.setfields('Profile_time')
    nrecs, _, _, _, _ = pid.inquire()
    time_since = pid.read(nRec=nrecs)
    pid.detach()
    time_since = np.asarray(time_since)

    # Get time of first data point
    time_seconds = time_since + time

    # Check whether in Greenland
    lat_mask = (latitude > roi.bounds.iloc[0][1]) & (latitude < roi.bounds.iloc[0][3])
    long_mask = (longitude < roi.bounds.iloc[0][2]) & (longitude > roi.bounds.iloc[0][0])

    valid = (lat_mask == True) & (long_mask == True)
    valid_lat, valid_lon = np.where(valid)

    valid_lon = longitude[valid]
    valid_lat = latitude[valid]
    valid_snow = snowfall_rate[valid]
    valid_time = time_seconds[valid]
    
    # Append
    lon_data.append(valid_lon)
    lat_data.append(valid_lat)

# Flatten lists
flat_lon = np.concatenate(lon_data).ravel()
flat_lat = np.concatenate(lat_data).ravel()

# Make a GeoDataFrame
gdf = gpd.GeoDataFrame(list(zip(flat_lon, flat_lat)))
gdf.rename(columns={0 : "lon", 1 : "lat", 2 : "snow", 3 : "time"}, inplace=True)
gdf['geometry'] = gpd.points_from_xy(gdf['lon'], gdf['lat'])
gdf.set_geometry(col='geometry', inplace=True)
gdf.crs = "EPSG:4326"
    
# Save
gdf.to_file(fp + 'cloudsat-greenland.shp')
    