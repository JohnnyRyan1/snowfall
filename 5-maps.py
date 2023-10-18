#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" 

DESCRIPTION

1. Maps

"""

# Import modules
import xarray as xr
import numpy as np
import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.colors import ListedColormap

# Define user
user = 'jryan4'

# Define path
path = '/Users/' + user + '/Dropbox (University of Oregon)/research/snowfall/data/'

# Define save path
savepath = '/Users/' + user + '/Dropbox (University of Oregon)/research/snowfall/manuscript/figures/'

# Import data
all_stats = pd.read_csv(path + 'all_stats.csv')

# Define years and regions
years = np.arange(2001, 2022)
regions = np.arange(1, 9, 1)

year=2010


#%%

# Import MERRA-2 grid
merra = xr.open_dataset(path + 'merra_resample/merra_' + str(year) + '.nc')

# Define empty array
sig_t = np.zeros((merra['t2m'].shape[1], merra['t2m'].shape[2]))
sig_s = np.zeros((merra['t2m'].shape[1], merra['t2m'].shape[2]))
sen_t = np.zeros((merra['t2m'].shape[1], merra['t2m'].shape[2]))
sen_s = np.zeros((merra['t2m'].shape[1], merra['t2m'].shape[2]))

#%%

# Assign values to grid
for r in range(all_stats.shape[0]):
    sig_t[all_stats['grid_cell_j'].iloc[r]][all_stats['grid_cell_i'].iloc[r]] = all_stats['p_values_t'].iloc[r]
    sig_s[all_stats['grid_cell_j'].iloc[r]][all_stats['grid_cell_i'].iloc[r]] = all_stats['p_values_s'].iloc[r]
    sen_t[all_stats['grid_cell_j'].iloc[r]][all_stats['grid_cell_i'].iloc[r]] = all_stats['coeffs_t'].iloc[r]
    sen_s[all_stats['grid_cell_j'].iloc[r]][all_stats['grid_cell_i'].iloc[r]] = all_stats['coeffs_s'].iloc[r]

sig_t[sig_t == 0] = np.nan
sig_s[sig_s == 0] = np.nan
sen_t[sen_t == 0] = np.nan
sen_s[sen_s == 0] = np.nan

#%%


# Save as NetCDF
ds_data = xr.Dataset(
data_vars={
    "sig_t": (("y", "x"), sig_t.astype('float32')),
    "sig_s": (("y", "x"), sig_s.astype('float32')),
    "sen_t": (("y", "x"), sen_t.astype('float32')),
    "sen_s": (("y", "x"), sen_s.astype('float32')),

},

coords={
    "y": (('y',), merra['y'].values),
    "x": (('x',), merra['x'].values),    
},

attrs={
    "Produced": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
    "Author":'Johnny Ryan', 
    "Email":'jryan4@uoregon.edu'
},
)

#%%
squares = sig_t
squares[np.isfinite(squares)] = 1
squares[np.isnan(squares)] = 0

#%%

# make a color map of fixed colors
cmap = colors.ListedColormap(['lightgrey', 'darkred'])
bounds=[0,0.5,1]
norm = colors.BoundaryNorm(bounds, cmap.N)

my_cmap = cmap(np.arange(cmap.N))

# Set alpha
my_cmap[:,-1] = np.linspace(0.4, 1, cmap.N)

# Create new colormap
my_cmap = ListedColormap(my_cmap)

# make a color map of fixed colors
cmap = colors.ListedColormap(['white', 'yellow'])
bounds=[0,0.5,1]
norm = colors.BoundaryNorm(bounds, cmap.N)

second_cmap = cmap(np.arange(cmap.N))

# Set alpha
second_cmap[:,-1] = np.linspace(0, 1, cmap.N)

# Create new colormap
second_cmap = ListedColormap(second_cmap)

# make a color map of fixed colors
cmap = colors.ListedColormap(['white', 'orange'])
bounds=[0,0.5,1]
norm = colors.BoundaryNorm(bounds, cmap.N)

third_cmap = cmap(np.arange(cmap.N))

# Set alpha
third_cmap[:,-1] = np.linspace(0, 1, cmap.N)

# Create new colormap
third_cmap = ListedColormap(third_cmap)

#%%

# Read stats DataFrame
all_stats = pd.read_csv(path + 'all_stats.csv')

# Identify three indexes in different regions
#ids = [82, 89, 126]
ids = [85, 90, 162]

# Define zero array
empty = np.zeros((merra['t2m'].shape[1], merra['t2m'].shape[2]))

# Assign value
empty[all_stats['grid_cell_j'].iloc[ids[0]], all_stats['grid_cell_i'].iloc[ids[0]]] = 1
empty[all_stats['grid_cell_j'].iloc[ids[1]], all_stats['grid_cell_i'].iloc[ids[1]]] = 1
empty[all_stats['grid_cell_j'].iloc[ids[2]], all_stats['grid_cell_i'].iloc[ids[2]]] = 1

#%%

# Identify three indexes in different regions
ids = [115, 21, 198, 101, 64, 20]

# N, N, NE, SW, CW, NW

# Define zero array
empty1 = np.zeros((merra['t2m'].shape[1], merra['t2m'].shape[2]))

# Assign value
empty1[all_stats['grid_cell_j'].iloc[ids[0]], all_stats['grid_cell_i'].iloc[ids[0]]] = 1
empty1[all_stats['grid_cell_j'].iloc[ids[1]], all_stats['grid_cell_i'].iloc[ids[1]]] = 1
empty1[all_stats['grid_cell_j'].iloc[ids[2]], all_stats['grid_cell_i'].iloc[ids[2]]] = 1
empty1[all_stats['grid_cell_j'].iloc[ids[3]], all_stats['grid_cell_i'].iloc[ids[3]]] = 1
empty1[all_stats['grid_cell_j'].iloc[ids[4]], all_stats['grid_cell_i'].iloc[ids[4]]] = 1
empty1[all_stats['grid_cell_j'].iloc[ids[5]], all_stats['grid_cell_i'].iloc[ids[5]]] = 1

#%%

fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(15,15), layout='constrained')
ax1.set_aspect('equal', adjustable='box')

ax1.pcolormesh(ds_data['x'], ds_data['y'], squares, edgecolor='k',
               cmap=my_cmap, norm=norm)
ax1.pcolormesh(ds_data['x'], ds_data['y'], empty, edgecolor='k',
               cmap=second_cmap, norm=norm)
ax1.pcolormesh(ds_data['x'], ds_data['y'], empty1, edgecolor='k',
               cmap=third_cmap, norm=norm)
ax1.tick_params(labelbottom=False)   
ax1.tick_params(labelleft=False)   
fig.savefig(savepath + 'fig_sx_maps.pdf')

#%%

###############################################################################
# Region version
###############################################################################

# Define empty array
sig_t = np.zeros((merra['t2m'].shape[1], merra['t2m'].shape[2]))

# Assign values to grid
for r in range(all_stats.shape[0]):
    sig_t[all_stats['grid_cell_j'].iloc[r]][all_stats['grid_cell_i'].iloc[r]] = all_stats['region'].iloc[r]

sig_t[sig_t == 0] = np.nan


#%%


# Save as NetCDF
ds_data = xr.Dataset(
data_vars={
    "sig_t": (("y", "x"), sig_t.astype('float32')),

},

coords={
    "y": (('y',), merra['y'].values),
    "x": (('x',), merra['x'].values),    
},

attrs={
    "Produced": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
    "Author":'Johnny Ryan', 
    "Email":'jryan4@uoregon.edu'
},
)

#%%
squares = sig_t
squares[np.isnan(squares)] = 0

fig, (ax1) = plt.subplots(nrows=1, ncols=1, figsize=(15,15), layout='constrained')
ax1.set_aspect('equal', adjustable='box')

ax1.pcolormesh(ds_data['x'], ds_data['y'], squares, edgecolor='k',
               cmap='gist_ncar_r')

ax1.tick_params(labelbottom=False)   
ax1.tick_params(labelleft=False)   
fig.savefig(savepath + 'fig_sx_region_maps.pdf')

























#%%