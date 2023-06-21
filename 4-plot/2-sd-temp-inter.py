#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" 

DESCRIPTION

1. Plot showing how snow depth, cumulative temperature, and timing of glacier 
ice exposure vary interannually.

2. Plot correlations between glacier ice exposure timing vs. snowdepth and June temp.

"""

#%%

# Import modules
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from scipy import stats
from matplotlib.offsetbox import AnchoredText
import statsmodels.api as sm
import pandas as pd

# Define path
path = '/Users/jryan4/Dropbox (University of Oregon)/research/snowfall/data/'

# Define save path
savepath = '/Users/jryan4/Dropbox (University of Oregon)/research/snowfall/manuscript/figures/'

# Import ISMIP 1 km grid
ismip_1km = xr.open_dataset(path + 'masks/1km-ISMIP6.nc')

# Define pixel of interest
x, y = 782, 522

exposures = []
coverages = []
pdds = []
sds = []
max_sd = []
max_sd_day = []

# Define years
years = np.arange(2001, 2022)

#%%

for year in years:
        
    print('Processing... %.0f' %year)
    
    # Import datasets
    modis = xr.open_dataset(path + 'modis_exposure_v2/MOD10A1_albedo_dates_' + str(year) + '.nc', engine='netcdf4')
    merra_sd = xr.open_dataset(path + 'merra_snowdepth/sd_' + str(year) + '.nc')
    merra_t2m = xr.open_dataset(path + 'merra_t2m/t2m_' + str(year) + '.nc') 
    
    # Produce glacier ice masks
    ice_mask_first_55 = modis['first_55'].values.astype(float)
    ice_mask_first_55[ismip_1km['ICE'] == 0] = np.nan
    ice_mask_albedo = modis['albedo'].values.astype(float)
    ice_mask_albedo[ismip_1km['ICE'] == 0] = np.nan
    ice_mask_elevation = ismip_1km['SRF'].values
    ice_mask_last_55 = modis['last_55'].values.astype(float)
    ice_mask_last_55[ismip_1km['ICE'] == 0] = np.nan
   
    # Find timing of exposure
    expo_test = ice_mask_first_55[x, y]
        
    # Find which MERRA grid cell contains this grid cell
    modis_point = [modis['latitude'].values[x, y],
                   modis['longitude'].values[x, y]]
    
    abslat = np.abs(merra_sd.latitude - modis_point[0])
    abslon = np.abs(merra_sd.longitude - modis_point[1])
    c = np.maximum(abslon, abslat)
    ([xloc], [yloc]) = np.where(c == np.min(c))
    
    # Convert temperature to positive degree days
    pdd = merra_t2m['t2m'][:, yloc, xloc] - 273
    pdd[pdd < 0] = 0
    
    # Compute timing of max snow depth
    max_sd.append(merra_sd['t2m'][:, yloc, xloc].values[merra_sd['t2m'][:, yloc, xloc].argmax().values])
    max_sd_day.append(merra_sd['time'].dt.dayofyear[[merra_sd['t2m'][:, yloc, xloc].argmax().values]].values[0])

    exposures.append(expo_test)
    coverages.append(ice_mask_last_55[x,y])
    pdds.append(np.array(pdd))
    sds.append(merra_sd['t2m'][:, yloc, xloc].values)

#%%

# Find some means
exposures = np.array(exposures)
exposures[exposures == 0] = np.nan
mean_expo = np.nanmedian(exposures)

coverages = np.array(coverages)
coverages[coverages == 0] = np.nan
mean_cover = np.nanmedian(coverages)

cum_pdds = np.cumsum(np.array(pdds), axis=1)
pdds_mean = np.mean(cum_pdds, axis=0)

sds = np.array(sds)
sds_mean = np.mean(sds, axis=0)

pdds = np.array(pdds)
mean_pdds_june = np.mean(pdds[:, 92:122], axis=1)

#%%
# Plot figure
fig, (ax1) = plt.subplots(nrows=1, ncols=1, figsize=(8, 4), layout='constrained')

ax1.plot(merra_sd['time'].dt.dayofyear, sds_mean, color='#0000a7', 
         lw=3, alpha=0.8, zorder=3)
ax2 = ax1.twinx()
ax2.plot(merra_t2m['time'].dt.dayofyear, pdds_mean, lw=3, color='#c1272d', 
         alpha=0.8, zorder=3)
ax1.axvline(x=mean_expo, ls='dashed', lw=2, color='k')
ax1.axvspan(np.nanpercentile(exposures, 10), np.nanpercentile(exposures, 90),
            color='grey', alpha=0.3)

# =============================================================================
# ax1.axvline(x=mean_cover, ls='dashed', lw=2, color='k')
# ax1.axvspan(np.nanpercentile(coverages, 10), np.nanpercentile(coverages, 90),
#             color='grey', alpha=0.3)
# =============================================================================

ax1.fill_between(merra_sd['time'].dt.dayofyear, np.percentile(sds, 10, axis=0), 
                 np.percentile(sds, 90, axis=0), 
                 alpha=0.2, color='#0000a7', zorder=2)

ax2.fill_between(merra_t2m['time'].dt.dayofyear, np.percentile(cum_pdds, 10, axis=0), 
                 np.percentile(cum_pdds, 90, axis=0),
                 alpha=0.2, color='#c1272d', zorder=2)

ax1.set_ylabel('Snow depth (m)', fontsize=14, color='#0000a7')
ax2.set_ylabel('Positive degree days (K)', fontsize=14, color='#c1272d')
ax1.set_xlabel('Day of year', fontsize=14)

ax1.set_ylim(0, 1.75)
ax2.set_ylim(0, 175)
ax1.set_xlim(100, 273)
ax1.grid(ls='dashed', lw=1, zorder=1)
ax1.tick_params(axis='both', which='major', labelsize=14)
ax2.tick_params(axis='both', which='major', labelsize=14)
ax1.tick_params(axis='y', colors='#0000a7')
ax2.tick_params(axis='y', colors='#c1272d')

ax1.text(0.38, 1.02, "Glacier ice exposure", fontsize=14, transform=ax1.transAxes)


fig.savefig(savepath + 'fig_sx_snow_depth_vs_temp_interannual.pdf')


#%%
mask = np.isfinite(exposures)

# Compute stats
e = exposures[mask]
p = mean_pdds_june[mask]
s = np.array(max_sd)[mask]

# Standardize
p = (p - np.mean(p)) / np.std(p)
s = (s - np.mean(s)) / np.std(s)

slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(p, e)
slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(s, e)

ts_slope1, ts_inter1, ts_low1, ts_high1 = stats.theilslopes(e, p, 0.90)
ts_slope2, ts_inter2, ts_low2, ts_high2 = stats.theilslopes(e, s, 0.90)

beta_pdd = r_value1 * (np.std(e) / np.std(p))
beta_sd = r_value2 * (np.std(e) / np.std(s))

# 
df = pd.DataFrame([list(p), list(s)]).T
df.columns = ['t','s']

X = df[['t', 's']]
y = pd.DataFrame(e, columns=['e'])

# Compute adjusted R2
r2_adj1 = 1 - (1 - r_value1**2) * (df.shape[0] - 1) / (df.shape[0] - df[['t']].shape[1] - 1)
r2_adj2 = 1 - (1 - r_value2**2) * (df.shape[0] - 1) / (df.shape[0] - df[['s']].shape[1] - 1)

# Multiple linear regression
X = sm.add_constant(X)
results = sm.OLS(y, X).fit()

def results_summary_to_df(results):

    pvals = results.pvalues
    coeff = results.params
    r2    = results.rsquared_adj

    results_df = pd.DataFrame({"pvals":pvals,
                               "coeff":coeff,
                               "r2_adj_both":r2,
                               "r2_temp":r2_adj1,
                               "r2_snow":r2_adj2,
                                })

    results_df = results_df[["coeff", "pvals", "r2_adj_both", "r2_temp", "r2_snow"]]
    return results_df

results_summary_to_df(results)

#%%

# Plot
fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(8, 4), sharey=True, layout='constrained')

ax1.scatter(p, e, color='#0000a7', zorder=2, s=75, alpha=0.6)
ax1.plot(p, p*slope1 + intercept1, color='k', lw=1)
#ax1.plot(p, (p*ts_slope1) + ts_inter1, color='r', lw=1)

ax2.scatter(s, e, color='#0000a7', zorder=2, s=75, alpha=0.6)
ax2.plot(s, s*slope2 + intercept2, color='k', lw=1)
#ax2.plot(s, s*ts_slope2 + ts_inter2, color='r', lw=1)

ax1.set_xlabel('Mean PDD June (K)', fontsize=14)
ax1.set_ylabel('Timing of exposure (DOY))', fontsize=14)
ax2.set_xlabel('Max snow depth (m)', fontsize=14)

ax1.tick_params(axis='both', which='major', labelsize=14)
ax2.tick_params(axis='both', which='major', labelsize=14)

ax1.grid(ls='dashed', lw=1, zorder=1)
ax2.grid(ls='dashed', lw=1, zorder=1)

# Add stats
textstr = '\n'.join((
    r'R$^{2}$ = %.2f' % (r_value1**2, ),))
text_box = AnchoredText(textstr, frameon=True, loc=3, pad=0.5, prop=dict(size=14))
text_box.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
plt.setp(text_box.patch, facecolor='white', alpha=0.7)
ax1.add_artist(text_box)

textstr = '\n'.join((
    r'R$^{2}$ = %.2f' % (r_value2**2, ),))
text_box = AnchoredText(textstr, frameon=True, loc=4, pad=0.5, prop=dict(size=14))
text_box.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
plt.setp(text_box.patch, facecolor='white', alpha=0.7)
ax2.add_artist(text_box)

for i, label in enumerate(list(years[mask])):
    ax1.annotate(label, (p[i], e[i]))
    ax2.annotate(label, (s[i], e[i]))

ax1.text(0.03, 0.89, "a", fontsize=24, transform=ax1.transAxes)
ax2.text(0.03, 0.89, "b", fontsize=24, transform=ax2.transAxes)

fig.savefig(savepath + 'fig_sx_sensitivity_analysis.pdf')


#%%

























