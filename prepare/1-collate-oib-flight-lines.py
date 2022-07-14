"""

Collate all OIB flight lines to single shapefile.

"""

# Import packages
import geopandas as gpd
import glob
import pathlib

# Define filepaths
fp = '/Users/jryan4/Dropbox (University of Oregon)/research/snowfall/data/oib-flight-lines/'
sp = '/Users/jryan4/Dropbox (University of Oregon)/research/snowfall/data/'

# Define file list
files = glob.glob(fp + '*.shp')

# Add to list
lines = []
name = []
for i in files:
    # Read file
    line = gpd.read_file(i)
    
    # Add to list
    lines.append(line['geometry'].values[0])
    name.append(pathlib.Path(i).stem)

# Make a new geodataframe
gdf = gpd.GeoDataFrame(list(zip(name, lines)))
gdf.rename(columns={0 : "name", 1 : "geometry"}, inplace=True)
gdf.set_geometry(col='geometry', inplace=True)
gdf.crs = "EPSG:4326"

# Save as file
gdf.to_file(sp + 'oib-flight-lines.shp')