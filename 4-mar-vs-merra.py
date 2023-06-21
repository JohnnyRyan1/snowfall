#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

MERRA-2 vs. MAR

"""

# Import modules
import xarray as xr
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error
from matplotlib.offsetbox import AnchoredText


# Define user
user = 'jryan4'

# Define path
path = '/Users/' + user + '/Dropbox (University of Oregon)/research/snowfall/data/'

# Define save path
savepath = '/Users/' + user + '/Dropbox (University of Oregon)/research/snowfall/manuscript/figures/'


# Read stats DataFrame
all_stats = pd.read_csv(path + 'all_stats.csv')

# Define years
years = np.arange(2001, 2022)

#%%

mar_t_list, mer_t_list = [], []
mar_s_list, mer_s_list = [], []

for year in years:
    
    print('Processing... %.0f' %year)
    
    # Read data
    mar = xr.open_dataset(path + 'mar_resample/mar_' + str(year) + '.nc')
    merra = xr.open_dataset(path + 'merra_resample/merra_' + str(year) + '.nc')
    
    mar_temp, merra_temp = [], []
    mar_snow, merra_snow = [], []
    
    # Loop over every MERRA grid cell
    for cell in range(all_stats.shape[0]):
        
        i = all_stats['grid_cell_i'].iloc[cell]
        j = all_stats['grid_cell_j'].iloc[cell]
            
        # Mean June air temperature
        merra_t = merra['t2m'][:, j, i] - 273
        mar_t = mar['t2m'][:, j, i]
          
        # Max snow depths
        merra_sd = merra['sd'][:, j, i].values[merra['sd'][:, j, i].argmax().values]
        
        if np.isfinite(mar['sd2'][:, j, i].values).sum() > 0:
            # Append
            mar_temp.append(np.mean(mar_t[151:181].values))
            merra_temp.append(np.mean(merra_t[92:122].values))
            
            mar_sd = mar['sd2'][:, j, i].values[mar['sd2'][:, j, i].argmax().values]
    
            # Append
            mar_snow.append(mar_sd)
            merra_snow.append(merra_sd)
            
        else:
            pass
    
    mar_t_list.append(mar_temp)
    mer_t_list.append(merra_temp)
    mar_s_list.append(mar_snow)
    mer_s_list.append(merra_snow)
    
#%%

# Convert to DataFrame
df_mar_t = pd.DataFrame(mar_t_list).T
df_mer_t = pd.DataFrame(mer_t_list).T

df_mar_s = pd.DataFrame(mar_s_list).T
df_mer_s = pd.DataFrame(mer_s_list).T

#%%

df_mar_t = df_mar_t[0:235]
df_mer_t = df_mer_t[0:235]
df_mar_s = df_mar_s[0:235]
df_mer_s = df_mer_s[0:235]

slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(df_mar_t.values.flatten(), df_mer_t.values.flatten())
slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(df_mar_s.values.flatten(), df_mer_s.values.flatten())

rms1 = mean_squared_error(df_mar_t.values.flatten(), df_mer_t.values.flatten(), squared=False)
rms2 = mean_squared_error(df_mar_s.values.flatten(), df_mer_s.values.flatten(), squared=False)

# Plot
fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(11, 4), 
                                      layout='constrained')
ax1.scatter(df_mar_t, df_mer_t, s=50, alpha=0.1, color='#0000a7', zorder=2)
ax2.scatter(df_mar_s, df_mer_s, s=50, alpha=0.1, color='#0000a7', zorder=2)

#ax1.hist2d(df_mar_t.values.flatten(), df_mer_t.values.flatten(), bins=30, zorder=2)
#ax2.hist2d(df_mar_s.values.flatten(), df_mer_s.values.flatten(), bins=30, zorder=2)

ax1.set_xlabel('Mean June air temperature (C) [MAR]', fontsize=14)
ax2.set_xlabel('Max. snow depth (m) [MAR]', fontsize=14)

ax1.set_ylabel('Mean June air temp. (C) [MERRA-2]', fontsize=14)
ax2.set_ylabel('Max. snow depth (m) [MERRA-2]', fontsize=14)

ax1.set_xlim(-14, 12)
ax1.set_ylim(-14, 12)
ax2.set_xlim(0, 20)
ax2.set_ylim(0, 12)
    
for ax in [ax1, ax2]:
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.grid(ls='dashed', lw=1, zorder=1)

# Add stats
textstr = '\n'.join((
    r'R$^{2}$ = %.2f' % (r_value1**2),
    r'RMSD = %.1f C' % rms1))
text_box = AnchoredText(textstr, frameon=True, loc=4, pad=0.5, prop=dict(size=13))
text_box.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
plt.setp(text_box.patch, facecolor='white', alpha=0.7)
ax1.add_artist(text_box)

textstr = '\n'.join((
    r'R$^{2}$ = %.2f' % (r_value2**2),
    r'RMSD = %.1f m' % rms2))
text_box = AnchoredText(textstr, frameon=True, loc=4, pad=0.5, prop=dict(size=13))
text_box.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
plt.setp(text_box.patch, facecolor='white', alpha=0.7)
ax2.add_artist(text_box)

ax1.text(0.01, 0.89, "a", fontsize=24, transform=ax1.transAxes)
ax2.text(0.01, 0.85, "b", fontsize=24, transform=ax2.transAxes)

fig.savefig(savepath + 'fig_sx_mar_vs_merra.pdf')


#%%
    
    
    
    
    
    
    
    
    
    


