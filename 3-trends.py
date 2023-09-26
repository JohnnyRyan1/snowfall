#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
DESCRIPTION

Determine whether or not there have been any trends in snow depth or air temperature

"""


#%%

# Import modules
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import norm, mstats
from matplotlib.offsetbox import AnchoredText


# Define user
user = 'jryan4'

# Define path
path = '/Users/' + user + '/Dropbox (University of Oregon)/research/snowfall/data/'

# Define save path
savepath = '/Users/' + user + '/Dropbox (University of Oregon)/research/snowfall/revision/'

# Read stats DataFrame
all_stats = pd.read_csv(path + 'all_stats.csv')


#%%
# Define years
years = np.arange(1981, 2022)

for year in years:
    
    print('Processing... %.0f' %year)
    
    # Import MERRA-2 grids
    merra_curr = xr.open_dataset(path + 'merra_resample/merra_' + str(year) + '.nc')
    merra_prev = xr.open_dataset(path + 'merra_resample/merra_' + str(year - 1) + '.nc')

    # Define empty lists
    spr_temp, jun_temp, win_temp = [], [], []
    spr_snow, snow_sum = [], []
    elevation, region = [], []
    
    # Loop over every MERRA grid cell
    for cell in range(all_stats.shape[0]):
        
        i = all_stats['grid_cell_i'].iloc[cell]
        j = all_stats['grid_cell_j'].iloc[cell]
        
        # Region and elevation
        region.append(all_stats['region'].iloc[cell])
        elevation.append(all_stats['elev'].iloc[cell])
        
        # Mean June air temperature
        t = merra_curr['t2m'][:, j, i] - 273
        jun_temp.append(np.mean(t[151:181]).values)
        
        # Mean Spring air temperature
        spr_temp.append(np.mean(t[59:151]).values)
        
        # Mean Spring snowfall
        spr_snow.append(merra_curr['sf'][:, j, i][59:151].sum().values)

        # Snowfall from Oct 1 to May 31
        sf_curr = merra_curr['sf'][:,j,i][0:151].sum().values
        sf_prev = merra_prev['sf'][:, j, i][273:].sum().values
        snow_sum.append(sf_curr + sf_prev)
        
        # Temperature from Oct 1 to May 31
        total = merra_curr['t2m'][:, j, i][0:151].sum().values + merra_prev['t2m'][:, j, i][273:].sum()
        win_temp.append((total/244).values)
        
    # Put in DataFrame
    df = pd.DataFrame(list(zip(region, elevation, 
                               spr_temp, jun_temp, win_temp, 
                               spr_snow, snow_sum)))
            
    df.columns = ['region', 'elevation', 
                  'spr_temp', 'jun_temp', 'win_temp', 
                  'spr_snow', 'snow_sum']
    
    # Combine with original index
    df['grid_cell_i'] = all_stats['grid_cell_i']
    df['grid_cell_j'] = all_stats['grid_cell_j']

    # Re-index
    df.set_index(['grid_cell_i', 'grid_cell_j'], inplace=True)
    
    # Save as csv
    df.to_csv(path + 'results/trends_' + str(year) + '.csv')

#%%

def mk_test(x, alpha = 0.05):  
    """   
    Input:
        x:   a vector of data
        alpha: significance level (0.05 default)

    Output:
        trend: tells the trend (increasing, decreasing or no trend)
        h: True (if trend is present) or False (if trend is absence)
        p: p value of the significance test
        z: normalized test statistics 

    Examples
    --------
      >>> x = np.random.rand(100)
      >>> trend,h,p,z = mk_test(x,0.05) 
    """
    n = len(x)

    # calculate S 
    s = 0
    for k in range(n-1):
        for j in range(k+1,n):
            s += np.sign(x[j] - x[k])

    # calculate the unique data
    unique_x = np.unique(x)
    g = len(unique_x)

    # calculate the var(s)
    if n == g: # there is no tie
        var_s = (n*(n-1)*(2*n+5))/18
    else: # there are some ties in data
        tp = np.zeros(unique_x.shape)
        for i in range(len(unique_x)):
            tp[i] = sum(unique_x[i] == x)
        var_s = (n*(n-1)*(2*n+5) + np.sum(tp*(tp-1)*(2*tp+5)))/18

    if s>0:
        z = (s - 1)/np.sqrt(var_s)
    elif s == 0:
            z = 0
    elif s<0:
        z = (s + 1)/np.sqrt(var_s)

    # calculate the p_value
    p = 2*(1-norm.cdf(abs(z))) # two tail test
    h = abs(z) > norm.ppf(1-alpha/2) 

    if (z<0) and h:
        trend = 'decreasing'
    elif (z>0) and h:
        trend = 'increasing'
    else:
        trend = 'no trend'

    return trend, h, p, z

#%%

years = np.arange(1981, 2022)

def compile_dfs(column): # Add functionalilty for preserving regions and elevation
        
    dfs = []
    size = []
    for year in years:
        
        # Import data
        df1 = pd.read_csv(path +'results/trends_' + str(year) + '.csv', 
                          index_col=['grid_cell_i', 'grid_cell_j'])

        # Append
        dfs.append(df1[column])
        size.append(df1.shape[0])
    
    # Merge
    df_merge = pd.concat(dfs, join='inner', axis=1)
    
    # Get largest DataFrame
    idx = years[np.argmax(np.array(size))]
    df2 = pd.read_csv(path +'results/trends_' + str(idx) + '.csv', 
                      index_col=['grid_cell_i', 'grid_cell_j'])
    
    # Merge with elevation
    df_merge2 = pd.concat([df_merge, df2[['elevation', 'region']]], join='inner', axis=1) 
    
    return df_merge, df_merge2[['elevation', 'region']]

#%%

all_s, elev = compile_dfs('snow_sum')
all_j, elev = compile_dfs('jun_temp')
all_spr, elev = compile_dfs('spr_temp')
all_win, elev = compile_dfs('win_temp')
spr_snow, elev = compile_dfs('spr_snow')

all_s.reset_index(inplace=True)
all_j.reset_index(inplace=True)
all_spr.reset_index(inplace=True)
spr_snow.reset_index(inplace=True)
all_win.reset_index(inplace=True)

# Remove columns
all_s.drop(columns = ['grid_cell_i', 'grid_cell_j'], inplace=True)
all_j.drop(columns = ['grid_cell_i', 'grid_cell_j'], inplace=True)
all_spr.drop(columns = ['grid_cell_i', 'grid_cell_j'], inplace=True)
all_win.drop(columns = ['grid_cell_i', 'grid_cell_j'], inplace=True)
spr_snow.drop(columns = ['grid_cell_i', 'grid_cell_j'], inplace=True)

#%%

# Whole ablation zone
trend,h,p,z = mk_test(all_j.mean(axis=0),0.05)

trend,h,p,z = mk_test(all_s.mean(axis=0),0.05)

trend,h,p,z = mk_test(all_win.mean(axis=0),0.05)

#%%

output_spr_temp = []
output_win_temp = []
output_jun_temp = []

output_win_snow = []
output_spr_snow = []

# Compute coefficients
for i in range(all_s.shape[0]):
    
    trend,h,p,z = mk_test(all_spr.iloc[i].values,0.05)
    output_spr_temp.append([trend, h, p, z])
    
    trend,h,p,z = mk_test(all_win.iloc[i].values,0.05)
    output_win_temp.append([trend, h, p, z])
    
    trend,h,p,z = mk_test(all_j.iloc[i].values,0.05)
    output_jun_temp.append([trend, h, p, z])
    
    trend,h,p,z = mk_test(all_s.iloc[i].values,0.05)
    output_win_snow.append([trend, h, p, z])
    
    trend,h,p,z = mk_test(spr_snow.iloc[i].values,0.05)
    output_spr_snow.append([trend, h, p, z])

jun_temp_df = pd.DataFrame(output_jun_temp, columns=['trend', 'h', 'p', 'z'])
spr_temp_df = pd.DataFrame(output_spr_temp, columns=['trend', 'h', 'p', 'z'])
win_temp_df = pd.DataFrame(output_win_temp, columns=['trend', 'h', 'p', 'z'])

win_snow_df = pd.DataFrame(output_win_snow, columns=['trend', 'h', 'p', 'z'])
spr_snow_df = pd.DataFrame(output_spr_snow, columns=['trend', 'h', 'p', 'z'])

# Add regions
jun_temp_df['region'] = elev['region'].values
spr_temp_df['region'] = elev['region'].values
win_temp_df['region'] = elev['region'].values

win_snow_df['region'] = elev['region'].values
spr_snow_df['region'] = elev['region'].values

# Stats
jun_temp_df[jun_temp_df['trend'] == 'increasing']
spr_temp_df[spr_temp_df['trend'] == 'increasing']
win_temp_df[win_temp_df['trend'] == 'increasing']
win_snow_df[win_snow_df['trend'] == 'decreasing']
spr_snow_df[spr_snow_df['trend'] == 'decreasing']

# Regional stats
counts = spr_temp_df.groupby(by='region')['trend'].count()
spr_increases = spr_temp_df[spr_temp_df['trend'] == 'increasing'].groupby(by='region')['trend'].count()
spr_increases / counts

jun_increases = jun_temp_df[jun_temp_df['trend'] == 'increasing'].groupby(by='region')['trend'].count()
jun_increases / counts

win_increases = win_temp_df[win_temp_df['trend'] == 'increasing'].groupby(by='region')['trend'].count()
win_increases / counts

win_snow_decreases = win_snow_df[win_snow_df['trend'] == 'decreasing'].groupby(by='region')['trend'].count()
win_snow_decreases / counts

win_snow_increases = win_snow_df[win_snow_df['trend'] == 'increasing'].groupby(by='region')['trend'].count()
win_snow_increases / counts

#%%

# Compute anomaly
spr_temp_norm = all_spr.sub(all_spr.mean(axis=1), axis=0)
jun_temp_norm = all_j.sub(all_j.mean(axis=1), axis=0)
win_temp_norm = all_win.sub(all_win.mean(axis=1), axis=0)

win_snow_norm = all_s.sub(all_s.mean(axis=1), axis=0)

#%%

# =============================================================================
# # Get increasing rows
# t_increase = spr_norm[spr_norm.index.isin(list(spr_df[spr_df['trend'] == 'increasing'].index))]
# t_neutral = spr_norm[spr_norm.index.isin(list(spr_df[spr_df['trend'] != 'increasing'].index))]
# s_decrease = sno_norm[sno_norm.index.isin(list(sno_df[sno_df['trend'] == 'decreasing'].index))]
# s_neutral = sno_norm[sno_norm.index.isin(list(sno_df[sno_df['trend'] != 'decreasing'].index))]
# 
# =============================================================================

# Stats
slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(np.arange(1981, 2022), jun_temp_norm.mean().values)
slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(np.arange(1981, 2022), win_snow_norm.mean().values)
slope3, intercept3, r_value3, p_value3, std_err3 = stats.linregress(np.arange(1981, 2022), win_temp_norm.mean().values)


#%%

# Plot
fig, (ax2, ax1, ax3) = plt.subplots(nrows=3, ncols=1, figsize=(12, 8), 
                                      layout='constrained', sharex=True)

ax1.plot(np.arange(1981, 2022), jun_temp_norm.mean().values, color='#c1272d', lw=2)
ax1.plot(np.arange(1981, 2022), np.arange(1981, 2022)*slope1 + intercept1, 
         ls='dashed', lw=2, color='k')
ax1.fill_between(np.arange(1981, 2022), np.percentile(jun_temp_norm, 10, axis=0), 
                 np.percentile(jun_temp_norm, 90, axis=0), color='grey', alpha=0.5)

ax2.plot(np.arange(1981, 2022), win_snow_norm.mean().values*365, color='#0000a7', lw=2)
ax2.plot(np.arange(1981, 2022), np.arange(1981, 2022)*slope2 + intercept2, 
         ls='dashed', lw=2, color='k')
ax2.fill_between(np.arange(1981, 2022), np.percentile(win_snow_norm*365, 10, axis=0), 
                 np.percentile(win_snow_norm*365, 90, axis=0), color='grey', alpha=0.5)

ax3.plot(np.arange(1981, 2022), win_temp_norm.mean().values, color='#e34a33', lw=2)
ax3.plot(np.arange(1981, 2022), np.arange(1981, 2022)*slope1 + intercept1, 
         ls='dashed', lw=2, color='k')
ax3.fill_between(np.arange(1981, 2022), np.percentile(win_temp_norm, 10, axis=0), 
                 np.percentile(win_temp_norm, 90, axis=0), color='grey', alpha=0.5)

ax1.set_ylim(-4, 4)
#ax2.set_ylim(-2,2)
ax2.set_xlim(1981, 2021)
ax1.set_ylabel('June air temp. \n anomaly (K)', fontsize=14)
ax2.set_ylabel('Cumulative snowfall \n [Oct 1 to May 31] (m)', fontsize=14)
ax3.set_ylabel('Winter air temp. anomaly \n [Oct 1 to May 31] (K)', fontsize=14)

ax1.tick_params(axis='both', which='major', labelsize=14)
ax2.tick_params(axis='both', which='major', labelsize=14)
ax3.tick_params(axis='both', which='major', labelsize=14)

ax1.grid(ls='dashed', lw=1, zorder=1)
ax2.grid(ls='dashed', lw=1, zorder=1)
ax3.grid(ls='dashed', lw=1, zorder=1)

ax1.text(0.01, 0.85, "b", fontsize=24, transform=ax1.transAxes)
ax2.text(0.01, 0.85, "a", fontsize=24, transform=ax2.transAxes)
ax3.text(0.01, 0.85, "c", fontsize=24, transform=ax3.transAxes)

# Add stats
textstr = '\n'.join((
    r'trend = +%.2f K yr$^{-1}$' % (slope1),
    r'p = %.2f' % p_value1))
text_box = AnchoredText(textstr, frameon=True, loc=4, pad=0.5, prop=dict(size=13))
text_box.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
plt.setp(text_box.patch, facecolor='white', alpha=0.7)
ax1.add_artist(text_box)

# Add stats
textstr = '\n'.join((
    r'trend = %.2f m yr$^{-1}$' % (slope2),
    r'p = %.2f' % p_value2))
text_box = AnchoredText(textstr, frameon=True, loc=4, pad=0.5, prop=dict(size=13))
text_box.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
plt.setp(text_box.patch, facecolor='white', alpha=0.7)
ax2.add_artist(text_box)

# Add stats
textstr = '\n'.join((
    r'trend = +%.2f K yr$^{-1}$' % (slope3),
    r'p = %.3f' % p_value3))
text_box = AnchoredText(textstr, frameon=True, loc=4, pad=0.5, prop=dict(size=13))
text_box.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
plt.setp(text_box.patch, facecolor='white', alpha=0.7)
ax3.add_artist(text_box)


fig.savefig(savepath + 'fig_3_trends.png', dpi=300)























    
    
    
        

