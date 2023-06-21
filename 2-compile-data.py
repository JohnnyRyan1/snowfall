#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" 

DESCRIPTION

1. Compile years 

2. Preliminary analysis

"""

# Import modules
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import statsmodels.api as sm
from matplotlib.offsetbox import AnchoredText

# Define user
user = 'johnnyryan'

# Define path
path = '/Users/' + user + '/Dropbox (University of Oregon)/research/snowfall/data/'

# Define save path
savepath = '/Users/' + user + '/Dropbox (University of Oregon)/research/snowfall/manuscript/figures/'

# Define years and regions
years = np.arange(2001, 2022)
regions = np.arange(1, 9, 1)

#%%

def compile_dfs(column):
        
    dfs = []
    size = []
    for year in years:
        
        # Import data
        df1 = pd.read_csv(path +'results/results_' + str(year) + '.csv', 
                          index_col=['grid_cell_i', 'grid_cell_j'])
        
        # Filter grid cells with few counts
        df1 = df1[df1['first_55_count'] >= 1]
            
        # Append
        dfs.append(df1[column])
        size.append(df1.shape[0])
    
    # Merge
    df_merge = pd.concat(dfs, join='inner', axis=1)
    
    # Get largest DataFrame
    idx = years[np.argmax(np.array(size))]
    df2 = pd.read_csv(path +'results/results_' + str(idx) + '.csv', 
                      index_col=['grid_cell_i', 'grid_cell_j'])
    
    # Merge with elevation
    df_merge2 = pd.concat([df_merge, df2[['elevation', 'region', 'lat', 'lon']]], join='inner', axis=1) 
    
    return df_merge, df_merge2[['elevation', 'region', 'lat', 'lon']]


#%%

all_s, elev = compile_dfs('snow_sum')
all_jun, elev = compile_dfs('jun_pdd')
all_jm, elev = compile_dfs('jun_temp')
all_spr, elev = compile_dfs('spr_temp')
all_spr_snow, elev = compile_dfs('spr_snow')
all_win, elev = compile_dfs('win_temp')
all_e, elev = compile_dfs('first_55_median')

# Remove months when June PDD is zero
all_s = all_s[all_jun.mean(axis=1) != 0]
all_spr_snow = all_spr_snow[all_jun.mean(axis=1) != 0]
all_win = all_win[all_jun.mean(axis=1) != 0]
all_jm = all_jm[all_jun.mean(axis=1) != 0]
all_spr = all_spr[all_jun.mean(axis=1) != 0]
all_e = all_e[all_jun.mean(axis=1) != 0]
elev = elev[all_jun.mean(axis=1) != 0]
all_jun = all_jun[all_jun.mean(axis=1) != 0]

#%%
coeffs_t = []
coeffs_s = []
r_values_t = []
r_values_s = []
p_values_t = []
p_values_s = []
std_t = []
std_s = []
coeffs_t_per_k = []

# Compute coefficients
for i in range(all_e.shape[0]):
    
    # Standardize    
    p = (all_jun.iloc[i].values - np.mean(all_jun.iloc[i].values)) / np.std(all_jun.iloc[i].values)
    s = (all_s.iloc[i].values - np.mean(all_s.iloc[i].values)) / np.std(all_s.iloc[i].values)
    e = all_e.iloc[i].values
    
    slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(p, e)
    slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(s, e)
    
    # Put back into DataFrame
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
    
    # Compute non-standardized linear regression
    slope3, intercept3, r_value3, p_value3, std_err3 = stats.linregress(all_jm.iloc[i].values, 
                                                                        e)
    # Append
    coeffs_t.append(results.params[1])
    coeffs_s.append(results.params[2])
    r_values_t.append(r_value1**2)
    r_values_s.append(r_value2**2)
    p_values_t.append(results.pvalues[1])
    p_values_s.append(results.pvalues[2])
    std_t.append(np.std(all_jun.iloc[i].values))
    std_s.append(np.std(all_s.iloc[i].values))
    coeffs_t_per_k.append(slope3)

#%%
all_stats = pd.DataFrame([elev['lon'].values, elev['lat'].values, elev['elevation'].values, elev['region'].values, 
                          std_t, std_s, coeffs_t, coeffs_s, r_values_t, r_values_s, 
                          p_values_t, p_values_s, coeffs_t_per_k]).T
all_stats.columns = ['x', 'y', 'elev', 'region', 'std_t', 'std_s', 'coeffs_t', 'coeffs_s', 
                     'r_values_t', 'r_values_s', 'p_values_t', 'p_values_s', 
                     'coeffs_t_per_k']

# Get grid cell index
all_s_idx = all_s.reset_index()

# Add to stats DataFrame
all_stats['grid_cell_i'] = all_s_idx['grid_cell_i']
all_stats['grid_cell_j'] = all_s_idx['grid_cell_j']

# Save
all_stats.to_csv(path + 'all_stats.csv', index=False)


#%%

###############################################################################
# Goodness-of-fit
###############################################################################

# Percentage of p < 0.05?
sig_t = np.sum((all_stats['p_values_t'] < 0.05) & (all_stats['coeffs_t'] < 0))
sig_s = np.sum((all_stats['p_values_s'] < 0.05) & (all_stats['coeffs_s'] > 0))

percent_sig_t = sig_t / np.sum(all_stats.shape[0])
percent_sig_s = sig_s / np.sum(all_stats.shape[0])

# R2 of significant grid cells
sig = all_stats[(all_stats['p_values_t'] < 0.05) & (all_stats['p_values_s'] < 0.05)]

std_r2_t = np.std(all_stats['r_values_t'])
std_r2_s = np.std(all_stats['r_values_s'])

mean_r2_t = all_stats['r_values_t'].mean()
mean_r2_s = all_stats['r_values_s'].mean()

mean_std_t = all_stats['std_t'].mean()
mean_std_s = all_stats['std_s'].mean()

mean_co_t = all_stats['coeffs_t'].mean()
mean_co_s = all_stats['coeffs_s'].mean()

# Number of grid cells where T is greater than SF
all_stats[np.abs(all_stats['coeffs_t']) > all_stats['coeffs_s']]['region']

#%%

# Plot three good r2 values
all_stats['combined'] = all_stats['r_values_t'] + all_stats['r_values_s']

# Identify three indexes in different regions
ids = [85, 90, 135]

# Standardize    
t1 = (all_jun.iloc[ids[0]].values - np.mean(all_jun.iloc[ids[0]].values)) / np.std(all_jun.iloc[ids[0]].values)
t2 = (all_jun.iloc[ids[1]].values - np.mean(all_jun.iloc[ids[1]].values)) / np.std(all_jun.iloc[ids[1]].values)
t3 = (all_jun.iloc[ids[2]].values - np.mean(all_jun.iloc[ids[2]].values)) / np.std(all_jun.iloc[ids[2]].values)

s1 = (all_s.iloc[ids[0]].values - np.mean(all_s.iloc[ids[0]].values)) / np.std(all_s.iloc[ids[0]].values)
s2 = (all_s.iloc[ids[1]].values - np.mean(all_s.iloc[ids[1]].values)) / np.std(all_s.iloc[ids[1]].values)
s3 = (all_s.iloc[ids[2]].values - np.mean(all_s.iloc[ids[2]].values)) / np.std(all_s.iloc[ids[2]].values)

e1 = all_e.iloc[ids[0]]
e2 = all_e.iloc[ids[1]]
e3 = all_e.iloc[ids[2]]

slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(s1, e1)
slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(s2, e2)
slope3, intercept3, r_value3, p_value3, std_err3 = stats.linregress(s3, e3)
slope4, intercept4, r_value4, p_value4, std_err4 = stats.linregress(t1, e1)
slope5, intercept5, r_value5, p_value5, std_err5 = stats.linregress(t2, e2)
slope6, intercept6, r_value6, p_value6, std_err6 = stats.linregress(t3, e3)

x = np.arange(-2, 4, 1)
    
# Plot figure
fig, ((ax1, ax2, ax3), 
      (ax4, ax5, ax6)) = plt.subplots(nrows=2, ncols=3, figsize=(12, 6), 
                                      layout='constrained', sharey=True)

ax1.scatter(s1, e1, color='#0000a7', zorder=2, s=75, alpha=0.6)
ax1.plot(x, x*slope1 + intercept1, color='k', lw=1)
#ax1.fill_between(x, ((x-conf1[1])*slope7 + intercept7), 
#                 ((x+conf1[1])*slope8 + intercept8), color='grey', alpha=0.4)

ax2.scatter(s2, e2, color='#0000a7', zorder=2, s=75, alpha=0.6)
ax2.plot(x, x*slope2 + intercept2, color='k', lw=1)
#ax2.fill_between(x, ((x+conf2[0])*slope2 + intercept2), 
#                 ((x+conf2[1])*slope2 + intercept2), color='grey', alpha=0.4)

ax3.scatter(s3, e3, color='#0000a7', zorder=2, s=75, alpha=0.6)
#ax3.plot(x, x*slope3 + intercept3, color='k', lw=1)
#ax3.fill_between(x, ((x+conf3[0])*slope3 + intercept3), 
#                 ((x+conf3[1])*slope3 + intercept3), color='grey', alpha=0.4)

ax4.scatter(t1, e1, color='#c1272d', zorder=2, s=75, alpha=0.6)
ax4.plot(x, x*slope4 + intercept4, color='k', lw=1)
#ax4.fill_between(x, ((x+conf4[0])*slope4 + intercept4), 
#                ((x+conf4[1])*slope4 + intercept4), color='grey', alpha=0.4)

ax5.scatter(t2, e2, color='#c1272d', zorder=2, s=75, alpha=0.6)
ax5.plot(x, x*slope5 + intercept5, color='k', lw=1)
#ax5.fill_between(x, ((x+conf5[0])*slope5 + intercept5), 
#                 ((x+conf5[1])*slope5 + intercept5), color='grey', alpha=0.4)

ax6.scatter(t3, e3, color='#c1272d', zorder=2, s=75, alpha=0.6)
ax6.plot(x, x*slope6 + intercept6, color='k', lw=1)
#ax6.fill_between(x, ((x+conf6[0])*slope6 + intercept6), 
#                 ((x+conf6[1])*slope6 + intercept6), color='grey', alpha=0.4)

ax4.set_xlabel('Standardized PDD June (K)', fontsize=14)
ax5.set_xlabel('Standardized PDD June (K)', fontsize=14)
ax6.set_xlabel('Standardized PDD June (K)', fontsize=14)
ax1.set_xlabel('Standardized snowfall (m)', fontsize=14)
ax2.set_xlabel('Standardized snowfall (m)', fontsize=14)
ax3.set_xlabel('Standardized snowfall (m)', fontsize=14)

ax1.set_ylabel('Timing of exposure (DOY)', fontsize=14)
ax4.set_ylabel('Timing of exposure (DOY)', fontsize=14)

for ax in [ax1, ax2, ax3, ax4, ax5, ax6]:
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.grid(ls='dashed', lw=1, zorder=1)
    ax.set_xlim(-2, 2.7)
    ax.set_ylim(150, 210)

# Add stats
textstr = '\n'.join((
    r'R$^{2}$ = %.2f' % (r_value1**2),
    r'slope = %.1f' %all_stats['coeffs_s'][ids[0]]))
text_box = AnchoredText(textstr, frameon=True, loc=1, pad=0.5, prop=dict(size=13))
text_box.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
plt.setp(text_box.patch, facecolor='white', alpha=0.7)
ax1.add_artist(text_box)

textstr = '\n'.join((
    r'R$^{2}$ = %.2f' % (r_value2**2),
    r'slope = %.1f' %all_stats['coeffs_s'][ids[1]]))
text_box = AnchoredText(textstr, frameon=True, loc=1, pad=0.5, prop=dict(size=13))
text_box.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
plt.setp(text_box.patch, facecolor='white', alpha=0.7)
ax2.add_artist(text_box)

textstr = '\n'.join((
    r'R$^{2}$ = %.2f' % (r_value3**2),
    r'slope = %.1f' %all_stats['coeffs_s'][ids[2]]))
text_box = AnchoredText(textstr, frameon=True, loc=1, pad=0.5, prop=dict(size=13))
text_box.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
plt.setp(text_box.patch, facecolor='white', alpha=0.7)
ax3.add_artist(text_box)

textstr = '\n'.join((
    r'R$^{2}$ = %.2f' % (r_value4**2),
    r'slope = %.1f' %all_stats['coeffs_t'][ids[0]]))
text_box = AnchoredText(textstr, frameon=True, loc=4, pad=0.5, prop=dict(size=13))
text_box.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
plt.setp(text_box.patch, facecolor='white', alpha=0.7)
ax4.add_artist(text_box)

textstr = '\n'.join((
    r'R$^{2}$ = %.2f' % (r_value5**2),
    r'slope = %.1f' %all_stats['coeffs_t'][ids[1]]))
text_box = AnchoredText(textstr, frameon=True, loc=4, pad=0.5, prop=dict(size=13))
text_box.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
plt.setp(text_box.patch, facecolor='white', alpha=0.7)
ax5.add_artist(text_box)

textstr = '\n'.join((
    r'R$^{2}$ = %.2f' % (r_value6**2),
    r'slope = %.1f' %all_stats['coeffs_t'][ids[2]]))
text_box = AnchoredText(textstr, frameon=True, loc=4, pad=0.5, prop=dict(size=13))
text_box.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
plt.setp(text_box.patch, facecolor='white', alpha=0.7)
ax6.add_artist(text_box)

ax1.text(0.03, 0.89, "a", fontsize=24, transform=ax1.transAxes)
ax2.text(0.03, 0.87, "b", fontsize=24, transform=ax2.transAxes)
ax3.text(0.03, 0.89, "c", fontsize=24, transform=ax3.transAxes)
ax4.text(0.03, 0.87, "d", fontsize=24, transform=ax4.transAxes)
ax5.text(0.03, 0.89, "e", fontsize=24, transform=ax5.transAxes)
ax6.text(0.03, 0.87, "f", fontsize=24, transform=ax6.transAxes)

ax1.text(0.37, 1.02, "Southwest", fontsize=16, transform=ax1.transAxes)
ax2.text(0.30, 1.02, "Central West", fontsize=16, transform=ax2.transAxes)
ax3.text(0.45, 1.02, "North", fontsize=16, transform=ax3.transAxes)

fig.savefig(savepath + 'fig_sx_sensitivity_analysis.pdf')


#%%

###############################################################################
# Sensitivity
###############################################################################

# Coefficients
mean_coeff_t = all_stats['coeffs_t'].mean()
mean_coeff_s = all_stats['coeffs_s'].mean()

std_coeffs_t = np.std(all_stats['coeffs_t'])
std_coeffs_s = np.std(all_stats['coeffs_s'])

# Summary table
is_summary = pd.DataFrame([percent_sig_t, percent_sig_s, mean_r2_t, mean_r2_s,
                           mean_coeff_t, mean_coeff_s, std_coeffs_t, std_coeffs_s,
                           mean_std_t, mean_std_s]).T
is_summary.columns = ['percent_sig_t', 'percent_sig_s', 'mean_r2_t', 'mean_r2_s',
                      'mean_coeff_t', 'mean_coeff_s', 'std_coeffs_t', 'std_coeffs_s',
                      'mean_std_t', 'mean_std_s']

is_summary.to_csv(path + 'summary_stats.csv', index=False)

#%%

# Plot histogram of sensitivities

bins = np.arange(all_stats['coeffs_t'].min(), all_stats['coeffs_s'].max(), 0.75)

# Plot figure
fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(12, 4), sharey=True, layout='constrained')

ax1.hist(np.abs(all_stats['coeffs_t']), bins=bins, color='#c1272d', 
         lw=1.5, alpha=0.6, zorder=3, edgecolor='k', label='Air temp.')
ax2.hist(all_stats['coeffs_s'], bins=bins, color='#0000a7', 
         lw=1.5, alpha=0.6, zorder=3, edgecolor='k', label='Snow depth')

ax1.axvline(x=np.mean(np.abs(all_stats['coeffs_t'])), lw=2, ls='dashed', c='k', zorder=4)
ax2.axvline(x=np.mean(all_stats['coeffs_s']), lw=2, ls='dashed', c='k', zorder=4)

ax1.set_ylabel('Frequency', fontsize=14)
ax1.set_xlabel('Glacier ice exposure sensitivity to \n June air temperature (days)', fontsize=14)
ax2.set_xlabel('Glacier ice exposure sensitivity to \n max. snow depth (days)', fontsize=14)

ax1.set_ylim(0, 35)
ax1.set_xlim(-5, 10)
ax1.grid(ls='dashed', lw=1, zorder=1)
ax1.tick_params(axis='both', which='major', labelsize=14)

ax2.set_ylim(0, 35)
ax2.set_xlim(-5, 10)
ax2.grid(ls='dashed', lw=1, zorder=1)
ax2.tick_params(axis='both', which='major', labelsize=14)

ax1.text(0.03, 0.89, "a", fontsize=24, transform=ax1.transAxes)
ax2.text(0.03, 0.87, "b", fontsize=24, transform=ax2.transAxes)

fig.savefig(savepath + 'fig_sx_snow_depth_vs_temp_histogram.pdf')

#%%
stats_is = []

for i in range(all_stats.shape[0]):
    
    # R2 and coefficients
    stats_is.append([all_jun.iloc[i].mean(),
                     all_s.iloc[i].mean(),
                     all_e.iloc[i].mean(),
                     all_jun.iloc[i].std(),
                     all_s.iloc[i].std(),
                     all_e.iloc[i].std(),
                     all_stats['r_values_t'].iloc[i], 
                     all_stats['r_values_s'].iloc[i],
                     np.abs(all_stats['coeffs_t'].iloc[i]),
                     np.abs(all_stats['coeffs_s'].iloc[i]),
                     all_stats['p_values_t'].iloc[i],
                     all_stats['p_values_s'].iloc[i],
                     all_stats['elev'].iloc[i],
                     all_stats['region'].iloc[i]])

is_df = pd.DataFrame(stats_is, columns=['mean_t', 'mean_s', 'mean_e', 
                                        'std_t', 'std_s', 'std_e', 
                                        'r2_t', 'r2_s', 
                                        'slope_t', 'slope_s',
                                        'p_values_t', 'p_values_s', 'elev',
                                        'region'])

is_df['weight'] = (np.abs(is_df['slope_t']) - is_df['slope_s']) / np.abs(is_df['slope_t'])
is_df_sig = is_df[(is_df['p_values_t'] < 0.05) & (is_df['p_values_s'] < 0.05)]

fractions = []
for r in np.arange(1, 9, 1):
    fractions.append([is_df[is_df['region'] == r].shape[0],
                      is_df_sig[is_df_sig['region'] == r].shape[0]])
                      
frac_df = pd.DataFrame(fractions, columns=['total', 'sig'])
frac_df['total_percent'] = frac_df['total'] / is_df.shape[0]
frac_df['sig_percent'] = frac_df['sig'] / is_df_sig.shape[0]
region_names = ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW']
frac_df['regions'] = region_names

#%%
###############################################################################
# Regional
###############################################################################

stats_region = []
for r in np.arange(1, 9, 1):
    
    # Filter region
    region_stats = all_stats[all_stats['region'] == r]
    
    t_region_stats = all_jun[elev['region'] == r]
    s_region_stats = all_s[elev['region'] == r]
    e_region_stats = all_e[elev['region'] == r]
    
    # R2 and coefficients
    stats_region.append([region_stats.shape[0],
                         region_stats['p_values_t'][region_stats['p_values_t'] < 0.05].shape[0]/region_stats.shape[0],
                         region_stats['p_values_s'][region_stats['p_values_s'] < 0.05].shape[0]/region_stats.shape[0],
                         e_region_stats.mean().mean(),
                         region_stats['r_values_t'].mean(), 
                         region_stats['r_values_s'].mean(),
                         np.abs(region_stats['coeffs_t'].mean()),
                         np.abs(region_stats['coeffs_s'].mean())])

region_names = ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW']
region_df = pd.DataFrame(stats_region, columns=['count', 'sig_t', 'sig_s', 'mean_e', 
                                                'r2_t', 'r2_s', 
                                                'slope_t', 'slope_s'])
region_df['region'] = region_names
#region_df['weight'] = (np.abs(region_df['slope_t']) - region_df['slope_s']) / np.abs(region_df['slope_t'])
region_df = region_df.T

region_df.to_csv(path + 'region_stats.csv', index=False)

#%%

slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(region_df['mean_s'], region_df['slope_t'])
slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(region_df['mean_s'], region_df['slope_s'])


# Some plots
fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(10, 4), 
                               layout='constrained')
ax1.scatter(region_df['mean_s'], region_df['slope_t'], s=100, alpha=0.6, color='#0000a7', zorder=2)
for i, txt in enumerate(region_names):
    ax1.annotate(txt, (region_df['mean_s'].iloc[i], region_df['slope_t'].iloc[i]), fontsize=14)

ax2.scatter(region_df['mean_s'], region_df['slope_s'], s=100, alpha=0.6, color='#0000a7', zorder=2)
for i, txt in enumerate(region_names):
    ax2.annotate(txt, (region_df['mean_s'].iloc[i], region_df['slope_s'].iloc[i]), fontsize=14)

ax1.set_xlabel('Mean max. snow depth 2001-2021 (m)', fontsize=14)
ax1.set_ylabel('Glacier ice exposure sensitivity to \n June air temperature (days)', fontsize=14)

ax2.set_xlabel('Mean max. snow depth 2001-2021 (m)', fontsize=14)
ax2.set_ylabel('Glacier ice exposure sensitivity to \n max. snow depth (days)', fontsize=14)
  
for ax in [ax1, ax2]:
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.grid(ls='dashed', lw=1, zorder=1)

ax1.set_ylim(1.8, 7.4)
ax2.set_ylim(0, 4.4)

# Add stats
textstr = '\n'.join((
    r'R$^{2}$ = %.2f' % (r_value1**2),
    r'p-value = %.3f' % p_value1))
text_box = AnchoredText(textstr, frameon=True, loc=1, pad=0.5, prop=dict(size=13))
text_box.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
plt.setp(text_box.patch, facecolor='white', alpha=0.7)
ax1.add_artist(text_box)

textstr = '\n'.join((
    r'R$^{2}$ = %.2f' % (r_value2**2),
    r'p-value = %.3f' % p_value2))
text_box = AnchoredText(textstr, frameon=True, loc=1, pad=0.5, prop=dict(size=13))
text_box.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
plt.setp(text_box.patch, facecolor='white', alpha=0.7)
ax2.add_artist(text_box)

ax1.text(0.01, 0.01, "a", fontsize=24, transform=ax1.transAxes)
ax2.text(0.01, 0.01, "b", fontsize=24, transform=ax2.transAxes)
    
fig.savefig(savepath + 'fig_sx_sensitivity_to_snow.pdf')



#%%

slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(region_df['std_t'], region_df['slope_s'])

fig, (ax1) = plt.subplots(nrows=1, ncols=1, figsize=(5, 4), layout='constrained')

ax1.scatter(region_df['mean_t'], region_df['slope_s'])
for i, txt in enumerate(region_names):
    ax1.annotate(txt, (region_df['mean_t'].iloc[i], region_df['slope_s'].iloc[i]), fontsize=14)



#%%
###############################################################################
# Relationship between 
#   1) Winter air temperature and winter snowfall
###############################################################################

t_vs_sf = []

# Compute coefficients
for i in range(all_win.shape[0]):
    
    # Define values    
    t = all_win.iloc[i].values
    s = all_s.iloc[i].values
    
    slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(t, s)
    
    t_vs_sf.append([p_value1, r_value1])

stats_df = pd.DataFrame(t_vs_sf, columns=['p', 'r'])
stats_df['region'] = all_stats['region']

# How many are significant during the study period?
stats_df[(stats_df['p'] < 0.05) & (stats_df['r'] > 0)]

sig = stats_df[(stats_df['p']< 0.05) & (stats_df['r'] > 0)]

percent = []
for r in regions:
    num_sig = sig[sig['region'] == r].shape[0]
    num_all = all_stats[all_stats['region'] == r].shape[0]
    fraction =  num_sig / num_all               
    percent.append([num_all, num_sig, fraction])

region_names = ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW']
t_vs_sf_df = pd.DataFrame(percent, columns=['num_sig', 'num_all', 'fraction'])
t_vs_sf_df['region'] = region_names


#%%

#%%
###############################################################################
# Relationship between 
#   1) June air temperature and cumulative snowfall
###############################################################################

t_vs_sf = []

# Compute coefficients
for i in range(all_jm.shape[0]):
    
    # Define values    
    t = all_jm.iloc[i].values
    s = all_s.iloc[i].values
    
    slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(t, s)
    
    t_vs_sf.append([p_value1, r_value1])

stats_df = pd.DataFrame(t_vs_sf, columns=['p', 'r'])
stats_df['region'] = all_stats['region']

# How many are significant during the study period?
stats_df[(stats_df['p'] < 0.05) & (stats_df['r'] > 0)]

sig = stats_df[(stats_df['p']< 0.05) & (stats_df['r'] > 0)]

percent = []
for r in regions:
    num_sig = sig[sig['region'] == r].shape[0]
    num_all = all_stats[all_stats['region'] == r].shape[0]
    fraction =  num_sig / num_all               
    percent.append([num_all, num_sig, fraction])

region_names = ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW']
t_vs_sf_df = pd.DataFrame(percent, columns=['num_sig', 'num_all', 'fraction'])
t_vs_sf_df['region'] = region_names

#%%















