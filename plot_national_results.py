#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_national_results.py
    python script to plot national HIA results from Shobha and Hai's datasets
Created on Mon Jul 11 10:42:48 2022
written by: Katelyn O'Dell
v0. initial code; 07.11.22
v1. final pre_submission
"""
#%% user inputs
# path to local project folder
prj_folder = '/Volumes/ODell_Files/Work_Files/GWU_072523/Ongoing_Projects/GeoXO/'

# if testing pct_reduce sensitivity, indicate here to load corect files
pct_reduce = '_0.30'#'_0.50' or '_0.15'

# path to datasets needed
# output from calc_annual_alert_hia.py
mort_data_path = '/Volumes/ODell_Files/from_pegasus/alerts_chronic_HIA_0.01_'
mort_version = 'v10' # v10 will be the final version once it finishes running
# output from calc_national_ervs.py
# add "_mort_" before all in the filename to load short-term mortality results reported in the supplement 
erv_data_path = prj_folder + 'HIA_results/erv_hia_results/national_alerts_HIA_'+pct_reduce[3:4]+'_mort_all_v7.nc'#'+ '_'+pct_reduce[3:]+'.nc'
# population output from prep_pop_data.py
pop_file = prj_folder + 'population_data/NASA_SEDAC_pop/regridded/regrided_2020pop_0.01_final.nc'
# airnow counts by census tract output from count_AQS_coverage.py
aqs_mon_path = prj_folder + 'concentration_data/EPA_AQI/AQS_obs_alert_counts_2020.csv'
# airnow counts by reporting area/ZCTAs output form count_airnowalert_coverage.py
airnow_path = prj_folder + 'misc/airnowalert_counts_dailymean.csv'

# path to datasets to plot/double check
# gridded baseline mortality rates output from grid_baseline_mort.py
gmort_br_path = '/Volumes/ODell_Files/from_pegasus/gridded_cnty_mort_final.csv'
# original baseline mortality rate values and shapefile
mort_br_path = prj_folder + 'health_data/baseline_rates/Underlying Cause of Death, 1999-2020.txt'
cnty_shp_fn = prj_folder + 'health_data/baseline_rates/tl_2019_us_county/tl_2019_us_county.shp'

# shapefiles for plotting
st_shp_fn = prj_folder + 'health_data/cb_2018_us_state_500k/cb_2018_us_state_500k.shp'
# shapefiles for plotting ZCTAs
zip_codes_shp_fn = prj_folder+'population_data/ACS_IPUMS/2020pop_2010zctas/nhgis0011_shape/nhgis0011_shapefile_tl2020_us_zcta_2010/US_zcta_2010_tl20.shp'

# local path to place figures
out_fig_path = prj_folder+'figures/ACX_PM/final/revisions/'

#%% import modules
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib as mplt
import geopandas as gpd
import pandas as pd

plt.style.use('seaborn-colorblind')

#%% load data
cnty_shp_df = gpd.read_file(cnty_shp_fn)
br_df = pd.read_csv(mort_br_path,sep='\t',nrows=3147)
gbr_df = pd.read_csv(gmort_br_path)

erv_fid = Dataset(erv_data_path)
geo_ervs = erv_fid['geo_erv'][:]
geo_erv_bm = erv_fid['geo_erv_bm'][:]
leo_erv_bm = erv_fid['abi1_erv_bm'][:]
erv_fid.close()

abi_fid = Dataset(mort_data_path + 'abi_og_'+mort_version+pct_reduce+'.nc')
abi_PM_popw = abi_fid['daily_PM_popw'][:]
abi_PMbm_popw = abi_fid['daily_PMbm_popw'][:]
abi_PM_areaw = abi_fid['daily_PM_areaw'][:]
abi_PMbm_areaw = abi_fid['daily_PMbm_areaw'][:]
abi_nobs = abi_fid['PM_obs'][:]
abi_nalerts = abi_fid['PM_alerts'][:]
abi_mort = abi_fid['base_mort'][:]
abi_mort_bm = abi_fid['bm_mort'][:]
abi_annavgPM = abi_fid['annavg_PM'][:]
abi_annavgPMbm = abi_fid['annavg_PM_bm'][:]
abi_area_covered = abi_fid['area_covered'][:]
abi_fid.close()

viirs_fid = Dataset(mort_data_path + 'viirs_og_'+mort_version+'_0.30.nc') # we use viirs og because there are times when viirs has data but abi does not
viirs_nobs = viirs_fid['PM_obs'][:]
viirs_nalerts = viirs_fid['PM_alerts'][:]
viirs_annavgPM = viirs_fid['annavg_PM'][:]
viirs_area_covered = viirs_fid['area_covered'][:]
viirs_mort = viirs_fid['base_mort'][:]
viirs_mort_bm = viirs_fid['bm_mort'][:]
viirs_annavgPM = viirs_fid['annavg_PM'][:]
viirs_annavgPMbm = viirs_fid['annavg_PM_bm'][:]
viirs_fid.close()

abi1pm_fid = Dataset(mort_data_path + 'abi1pm_proxy_'+mort_version+pct_reduce+'.nc') # change proxy to og to run the sensitivity analysis in Table S2
abi1pm_nobs = abi1pm_fid['PM_obs'][:]
abi1pm_nalerts = abi1pm_fid['PM_alerts'][:]
abi1pm_mort_bm = abi1pm_fid['bm_mort'][:]
abi1pm_annavgPM = abi1pm_fid['annavg_PM'][:]
abi1pm_mort = abi_fid['base_mort'][:]
abi1pm_annavgPMbm = abi1pm_fid['annavg_PM_bm'][:]
abi1pm_area_covered = abi1pm_fid['area_covered'][:]
abi1pm_fid.close()

an_df = pd.read_csv(airnow_path)
aqs_mon1 = pd.read_csv(aqs_mon_path)

nc_fid = Dataset(pop_file)
pop = nc_fid['population'][:].data
grid_area = nc_fid['grid_area'][:].data
area_mask = nc_fid['us_mask'][:].data
cal_area_mask = nc_fid['cal_mask'][:].data
lat = nc_fid['lat'][:].data
lon = nc_fid['lon'][:].data
nc_fid.close()
plt_mask = np.where(area_mask==0,np.nan,area_mask)

zipdf = gpd.read_file(zip_codes_shp_fn)

glon, glat = np.meshgrid(lon,lat)
day_nums = np.arange(0,366)

#%% initial figures/data checking
### FIGURES A and B: gridded baseline rate, compare to original file ###
# note, there will be a few slight differences because we replaced unreliable
# and suppressed counties in the original file with state-level baseline mort rates

# remove outside contig us fips from aqs_mon
inds_drop = []
aqs_mon1.reset_index(inplace=True,drop=True)
for i in range(aqs_mon1.shape[0]):
    fip = aqs_mon1['STATEFP10'].iloc[i]
    if fip in [2,15]:
        inds_drop.append(i)
aqs_mon = aqs_mon1.drop(inds_drop,axis=0)
aqs_mon.reset_index(inplace=True,drop=True)

# first adjust mismatched county code
ind_up = np.where(br_df['County Code'] == 46113)[0]
br_df.loc[ind_up,'County Code'] = 46102
# reshape gridded values
gmort = gbr_df['cnty_mort_rate'].values
pmort = gmort.reshape([2600, 5900])
# plot gridded values
'''
ax = plt.axes(projection=ccrs.PlateCarree())
plt.pcolormesh(glon, glat, pmort,vmin=0,vmax=1500)
ax.coastlines()
ax.set_title('gridded baseline mort')
#plt.show()

# now for original file to compare
# drop values outside the contig US
codes_drop = ['02','15','60','66','69','72','78']
for code_drop in codes_drop:
    cnty_shp_df.drop(cnty_shp_df[cnty_shp_df['STATEFP']==code_drop].index,inplace=True)
cnty_shp_df.reset_index(inplace=True,drop=True)

# add baseline mort rate to the shapefile
cnty_mort_rate = []
for ci in range(cnty_shp_df.shape[0]):
    cnty_id = int(cnty_shp_df['GEOID'].iloc[ci])
    ind = np.where(br_df['County Code']==cnty_id)[0]
    if len(ind)<1:
        # print(cnty_id, cnty_shp_df['NAMELSAD'].iloc[ci])
        cnty_mort_rate.append(-9999.9)
        continue
    br = br_df['Crude Rate'].iloc[ind[0]]
    if br in ['Unreliable','Suppressed']:
        cnty_mort_rate.append(-8888.8)
        continue
    cnty_mort_rate.append(float(br_df['Crude Rate'].iloc[ind[0]]))
cnty_shp_df['Crude Mort Rate'] = cnty_mort_rate
ax = cnty_shp_df.plot("Crude Mort Rate", legend=True,vmin=0,vmax=1500)
ax.set_axis_off()
ax.set_title('original county level baseline mort')

### FIGURE C: data in asthma ERVs file
fig = plt.figure()
ax = plt.axes(projection=ccrs.PlateCarree())
plt.pcolormesh(glon, glat, (10**5)*geo_ervs[1,:,:]/pop,vmin=0,vmax=10)
ax.coastlines()
ax.set_title('geo_ervs')
#plt.show()
'''
#%% FIGURE 1 - observations and annual averages
# prep data for plotting
# edit abi datasets so that values with no data on the map show up as gray for the annual means
abi1pm_annavg_plot = np.where(abi1pm_nobs==0,-99999,abi1pm_annavgPM)
abi_annavg_plot = np.where(abi_nobs==0,-99999,abi_annavgPM)
# merge zipcode shapefile with data file for AN-RA dataset 
zipdf = zipdf.astype({"GEOID10": np.int64})
zipoutdf = zipdf.merge(an_df, left_on='GEOID10', right_on='geoid10')
zipoutdf['obs_count_plot'] = np.where(pd.isna(zipoutdf['obs_count']),0,zipoutdf['obs_count'])
zipoutdf['pm_mean_plot'] = np.where(pd.isna(zipoutdf['pm_mean']),0,zipoutdf['pm_mean'])
# load state df and drop states outside us
st_df = gpd.read_file(st_shp_fn)
for state in ['AK','HI','AS','PR','VI','GU','MP']:
    st_df.drop(st_df[st_df['STUSPS'] ==state].index, inplace = True)
# allign projections for plotting
zipoutdf_plot = zipoutdf.to_crs(ccrs.PlateCarree())
st_df_plot = st_df.to_crs(ccrs.PlateCarree())
# create colormaps we want to use
my_magma = mplt.cm.get_cmap('plasma').copy()
my_magma.set_under('darkgray')
my_red = mplt.cm.get_cmap('Reds').copy()
my_red.set_under('darkgray')
# put titles in a list so we can add them to panels in a loop
titles = ['(a) AirNow-RA','(e) AirNow-RA','(b) VIIRS','(f) VIIRS',
          '(c) ABI-1pm', '(g) ABI-1pm','(d) ABI-Daytime','(h) ABI-Daytime']

# now, finally, create figure and plot these data
fig, axarr1 = plt.subplots(nrows=4,ncols=2,subplot_kw={'projection': ccrs.PlateCarree()},
                          figsize=(8,10))
axarr = axarr1.flatten()
# set up airnow panels and plot
for i in range(2):
    ax = st_df_plot.plot(edgecolor='gray',color='darkgray',ax = axarr[i],zorder=1)
    ax.axis("off")
    ax.set_title(titles[i],fontsize=18)
ax = zipoutdf_plot.plot('obs_count_plot', ax = axarr[0], vmin=1,vmax=366,cmap=my_magma, 
                   edgecolor='none',aspect='auto')
ax = zipoutdf_plot.plot('pm_mean_plot', ax = axarr[1], vmin=1,vmax=30,cmap=my_red, 
                   edgecolor='none',aspect='auto')

# set up satellite pannels and plot
for i in range(2,len(axarr)):
    axarr[i].axis("off")
    axarr[i].add_feature(cfeature.STATES.with_scale('50m'),edgecolor='gray')
    axarr[i].set_title(titles[i],fontsize=18)
cs1=axarr[2].pcolormesh(glon,glat,viirs_nobs*plt_mask,vmin=1,vmax=366,
            transform=ccrs.PlateCarree(),cmap=my_magma)
cs2=axarr[3].pcolormesh(glon,glat,viirs_annavgPM*plt_mask,vmin=1,vmax=30,
            transform=ccrs.PlateCarree(),cmap=my_red)
axarr[4].pcolormesh(glon,glat,abi1pm_nobs*plt_mask,vmin=1,vmax=366,
            transform=ccrs.PlateCarree(),cmap=my_magma)
axarr[5].pcolormesh(glon,glat,abi1pm_annavg_plot*plt_mask,vmin=1,vmax=30,
            transform=ccrs.PlateCarree(),cmap=my_red)
axarr[6].pcolormesh(glon,glat, abi_nobs*plt_mask,vmin=1,vmax=366,
            transform=ccrs.PlateCarree(),cmap=my_magma)
axarr[7].pcolormesh(glon,glat,abi_annavg_plot*plt_mask,vmin=1,vmax=30,
            transform=ccrs.PlateCarree(),cmap=my_red)
# make colorbars
cbar = fig.colorbar(cs1,ax=[axarr[-2]],orientation='horizontal',pad=0.1,shrink=0.6,ticks = [1,120,240,366])
cbar.set_label(label='Total Days of\nObservations [days]',size=16)
cbar = fig.colorbar(cs2,ax=[axarr[-1]],orientation='horizontal',pad=0.1,shrink=0.6,ticks = [1,10,20,30],extend='max')
cbar.set_label(label='Mean PM$_{2.5}$ [$\mu$g m$^{-3}$]',size=16)
# savefigure
plt.savefig(out_fig_path+'annual2020_observations_nat_all.png',dpi=300)
#plt.show()
plt.close()

# print numbers for figure 1: calculate average area covered by each dataset each day
mean_area_viirs = np.nansum((viirs_nobs/366.0)*grid_area*plt_mask)
mean_area_abiD = np.nansum((abi_nobs/366.0)*grid_area*plt_mask)
mean_area_abi1pm = np.nansum((abi1pm_nobs/366.0)*grid_area*plt_mask)
tot_grid_area = np.nansum(grid_area*plt_mask)

mean_area_an_zip = np.nansum((zipoutdf['ALAND10'].values + zipoutdf['AWATER10'].values)*zipoutdf['obs_count'].values)/(366.0*10**6)
tot_zip_area = np.nansum((zipoutdf['ALAND10'].values + zipoutdf['AWATER10'].values))/(1.0*10**6)

# annual mean PM
print('abi-Daytime',np.nanmean(abi_annavgPM*plt_mask))
print('abi-1pm',np.nanmean(abi1pm_annavgPM*plt_mask))
print('viirs',np.nanmean(viirs_annavgPM*plt_mask))
print('an',np.nanmean(zipoutdf['pm_mean']))

#%% Figure 2 - alert days maps
# make list of datasets to plot in a loop
data = [abi_nalerts*plt_mask,
        (viirs_nalerts)*plt_mask,
        (abi1pm_nalerts)*plt_mask]
# list of titles
titles = ['(a) ABI-Daytime','(b) VIIRS','(c) ABI-1pm','(d) AirNow-RA']
# create figure
fig, axarr = plt.subplots(nrows=2,ncols=2,subplot_kw={'projection': ccrs.PlateCarree()},
                          figsize=(11,5))
# loop through and plot satellite datasets
axarr=axarr.flatten()
for i in range(len(axarr[:-1])):
    axarr[i].axis("off")
    axarr[i].add_feature(cfeature.STATES.with_scale('50m'),edgecolor='gray')
    cs = axarr[i].pcolormesh(glon,glat,data[i],vmin=1,vmax=30,transform=ccrs.PlateCarree(),cmap=my_magma)
    axarr[i].set_title(titles[i],fontsize=12)
# plot airnow dataset
ax = st_df_plot.plot(edgecolor='gray',color='darkgray',ax = axarr[-1])
ax = zipoutdf_plot.plot('alerts_max', ax = axarr[-1], vmin=1,vmax=30,cmap=my_magma, 
                   edgecolor='none',aspect='auto', missing_kwds={'color': 'darkgrey'})
ax.axis("off")
ax.set_title(titles[-1],fontsize=12)
# adjust figure dimensions and add colobar  
plt.subplots_adjust(bottom=0.17)
cax = plt.axes([0.3, 0.1, 0.4, 0.04])
cbar = fig.colorbar(cs,cax=cax,orientation='horizontal',pad=0.1,shrink=0.3,
                ticks = [1,5,10,15,20,25,30],extend='max')
cbar.set_label(label='N Alert Days, 2020',size=12)
# save figure
plt.savefig(out_fig_path+'annual2020_alerts_nat_all.png',dpi=300)
#plt.show()
plt.close()


# Difference Plot for SI
titles = ['(a) ABI, Daytime - VIIRS','(b) ABI, Daytime - ABI, 1pm']
# create figure
fig, axarr = plt.subplots(nrows=2,ncols=1,subplot_kw={'projection': ccrs.PlateCarree()},
                          figsize=(6,5))
data = [(abi_nalerts-viirs_nalerts)*plt_mask,
        (abi_nalerts-abi1pm_nalerts)*plt_mask]
# loop through and plot satellite datasets
axarr=axarr.flatten()
for i in range(len(axarr)):
    axarr[i].axis("off")
    axarr[i].add_feature(cfeature.STATES.with_scale('50m'),edgecolor='gray')
    cs = axarr[i].pcolormesh(glon,glat,data[i],vmin = -10,vmax=10,
                             transform=ccrs.PlateCarree(),cmap='RdBu_r')
    axarr[i].set_title(titles[i],fontsize=12)
# plot airnow dataset
# adjust figure dimensions and add colobar  
plt.subplots_adjust(bottom=0.17)
cax = plt.axes([0.3, 0.1, 0.4, 0.04])
cbar = fig.colorbar(cs,cax=cax,orientation='horizontal',pad=0.1,shrink=0.3,
                ticks = [-10,-5,0,5,10],extend='both')
cbar.set_label(label='Difference in Alert Days, 2020',size=12)
# save figure
plt.savefig(out_fig_path+'annual2020_alerts_diff_all.png',dpi=300)
#plt.show()
plt.close()


# calculate and print numbers for figure discussion in the manuscript
print('mean alert days for each dataset')
print('ABI-Daytime',np.nanmean(abi_nalerts*plt_mask))
print('ABI-1pm',np.nanmean(abi1pm_nalerts*plt_mask))
print('VIIRS',np.nanmean(viirs_nalerts*plt_mask))

print('total alert days in each dataset')
print('ABI-Daytime',np.nansum(abi_nalerts*plt_mask))
print('ABI-1pm',np.nansum(abi1pm_nalerts*plt_mask))
print('VIIRS',np.nansum(viirs_nalerts*plt_mask))

# for zipcodes do an area weighted mean since zctas can be different sizes
print('AN alerts',np.nansum(zipoutdf_plot['alerts_max']*(zipoutdf_plot['AWATER10']+zipoutdf_plot['ALAND10']))/np.nansum(zipoutdf_plot['AWATER10']+zipoutdf_plot['ALAND10']))
print('AN alerts',np.nansum(zipoutdf['alerts_max']*(zipoutdf['pop']))/np.nansum(zipoutdf['pop']))

#%% FIGURE 3 - person alerts
# prepare data for calculating person alers
# first, make a plot mask with zeros isntead of nans for calculating pop-weighted means
us_mask = np.where(np.isnan(plt_mask),0,1)
pop_us = pop*us_mask
# calculate pop-weighted person-days in each dataset
# here an = AirNow-RA and mon = AirNow-MT 
abi_person_days = np.nansum(abi_nobs*pop_us)/np.nansum(pop_us)
abi1pm_person_days = np.nansum(abi1pm_nobs*pop_us)/np.nansum(pop_us)
viirs_person_days = np.nansum(viirs_nobs*pop_us)/np.nansum(pop_us)
an_person_days = np.nansum(an_df['obs_count'].values*an_df['pop'].values)/np.nansum(an_df['pop'].values)
mon_person_days = np.nansum(aqs_mon['obs_count'].values*aqs_mon['CL8AA2020'].values)/np.nansum(aqs_mon['CL8AA2020'].values)

# calc popw alerts
abi_popw_alerts = np.nansum(abi_nalerts*pop_us)/np.nansum(pop_us)
abi1pm_popw_alerts = np.nansum(abi1pm_nalerts*pop_us)/np.nansum(pop_us)
viirs_popw_alerts = np.nansum(viirs_nalerts*pop_us)/np.nansum(pop_us)
an_popw_alerts = np.nansum(an_df['alerts_max'].values*an_df['pop'].values)/np.nansum(an_df['pop'].values)

# calculate total person-alerts in each dataset
abi_person_alerts = np.nansum(abi_nalerts*pop_us)
abi1pm_person_alerts = np.nansum(abi1pm_nalerts*pop_us)
viirs_person_alerts = np.nansum(viirs_nalerts*pop_us)
an_person_alerts = np.nansum(an_df['alerts_max'].values*an_df['pop'].values)
mon_person_alerts = np.nansum(aqs_mon['CL8AA2020'].values*aqs_mon['alert_count'].values)
# define colors for each dataset
# these were picked using colorbrewer2: https://colorbrewer2.org/
colors = [#'#fdbf6f',# monitors
          '#ff7f00', # an
          '#b2df8a', # viirs
          '#a6cee3', # abi1pm
          '#1f78b4'] #abi
# put lables, person days and person alerts counts in lists for plotting
data_lables = ['AirNow-RA','VIIRS','ABI-1pm','ABI-Daytime']
person_days = [an_person_days, 
               viirs_person_days, abi1pm_person_days, abi_person_days]
person_alerts = [an_person_alerts,
                 viirs_person_alerts,abi1pm_person_alerts,abi_person_alerts]
# create figure and plot data
fig, ax = plt.subplots(1,2,figsize=(9,5))
ax[0].bar(data_lables, person_days, color = colors)
ax[1].bar(data_lables, person_alerts, color = colors)
ax[0].set_ylabel('Population-weighted\nN days with observations')
ax[1].set_ylabel('Person-Alerts')
ax[0].set_xticks(data_lables, data_lables,rotation=45)
ax[1].set_xticks(data_lables, data_lables,rotation=45)
plt.subplots_adjust(bottom = 0.3)
plt.savefig(out_fig_path + 'person_alerts_all.png',dpi=300)
plt.close()

# print numbers used in discussion of this figure
print('person-alert totals')
print(data_lables, np.round(person_alerts,-3))
print('popw mean days with obs')
print(data_lables,'\n',np.round(person_days,2))

# days in abi_daytime that are not in abi1pm


#%% FIGURE 4 - PM annual average diff
# create color maps
my_blues = mplt.cm.get_cmap('Blues')
my_blues.set_under('darkgray')
my_viridis = mplt.cm.get_cmap('viridis')
my_viridis.set_under('darkgray')

pre_nan_inds = np.where(np.isnan(abi_annavgPM))
# put data, colormaps, and titles in a list for plotting
diff1 = abi_annavgPM-abi_annavgPMbm
#diff2 = abi1pm_annavgPM-abi1pm_annavgPMbm
diff2 = abi_annavgPM-abi1pm_annavgPMbm
diff3 = diff1-diff2
diff1[pre_nan_inds] = -999 # set to -999 so it shows up as gray on the plot
diff2[pre_nan_inds] = -999 
diff3[pre_nan_inds] = -999 

data =  [diff1*plt_mask,diff2*plt_mask,diff3*plt_mask]
cmaps = [my_blues,my_blues,my_viridis]
titles = ['(a) ABI-Daytime','(b) ABI-1pm','(c) Difference']

# create figure and plot datasets
fig, axarr = plt.subplots(nrows=1,ncols=3,subplot_kw={'projection': ccrs.PlateCarree()},figsize=(11,3))
i = 0
for ax in axarr.flatten():
    ax.patch.set_visible(False)
    ax.axis("off")
    ax.add_feature(cfeature.STATES.with_scale('50m'),edgecolor='gray')
    cs = ax.pcolormesh(glon,glat,data[i],vmin=0,vmax=3,transform=ccrs.PlateCarree(),cmap=my_blues)
    cbar = fig.colorbar(cs,ax=ax,orientation='horizontal',pad=0,shrink=0.6,
                        ticks=[0.1,1,2,3],extend='max')
    cbar.set_label(label='PM$_{2.5}$ Reduction '+pct_reduce[1:],size=16) #+pct_reduce[1:]
    ax.set_title(titles[i],fontsize=18)
    i += 1
plt.tight_layout()
plt.savefig(out_fig_path+'annual_mean_pm25_reduction'+pct_reduce+'.png',dpi=350)
plt.show()
#plt.close()

# print numbers for figure discussion
# national mean pm2.5 exposure reductions
# remove -999s
diff1[pre_nan_inds] = np.nan # set to -999 so it shows up as gray on the plot
diff2[pre_nan_inds] = np.nan
diff3[pre_nan_inds] = np.nan
data =  [diff1*plt_mask,diff2*plt_mask,diff3*plt_mask]
print(np.nanmean(data[0]),np.nanmean(data[1]),np.nanmean(data[2]))
# population-weighted mean exposure reduction
pop_mean_base = np.nansum(pop_us*abi_annavgPM)/np.nansum(pop_us)
pop_mean_abidaytime = np.nansum(pop_us*abi_annavgPMbm)/np.nansum(pop_us)
pop_mean_abi1pm = np.nansum(pop_us*abi1pm_annavgPMbm)/np.nansum(pop_us)

# california specific change
cal_pop_mean_base = np.nansum(pop_us*abi_annavgPM*cal_area_mask)/np.nansum(pop_us*cal_area_mask)
cal_pop_mean_abidaytime = np.nansum(pop_us*abi_annavgPMbm*cal_area_mask)/np.nansum(pop_us*cal_area_mask)
cal_pop_mean_abi1pm = np.nansum(pop_us*abi1pm_annavgPMbm*cal_area_mask)/np.nansum(pop_us*cal_area_mask)

#%% FIGURE 5 - bar charts
# prepare data for plotting
# calcualte total, lower bound, and upper bound of ABI-Daytime reduced mortalities
tot_plot = np.nansum(abi_mort[0,:,:]*area_mask)-np.nansum(abi_mort_bm[0,:,:]*area_mask)
low_ci = np.nansum(abi_mort[1,:,:]*area_mask)-np.nansum(abi_mort_bm[1,:,:]*area_mask)
hi_ci = np.nansum(abi_mort[2,:,:]*area_mask)-np.nansum(abi_mort_bm[2,:,:]*area_mask)
# calculate the same for ABI-1pm
tot_plot1 = np.nansum(abi_mort[0,:,:]*area_mask)-np.nansum(abi1pm_mort_bm[0,:,:]*area_mask)
low_ci1 = np.nansum(abi_mort[1,:,:]*area_mask)-np.nansum(abi1pm_mort_bm[1,:,:]*area_mask)
hi_ci1 = np.nansum(abi_mort[2,:,:]*area_mask)-np.nansum(abi1pm_mort_bm[2,:,:]*area_mask)
# define function for converting deaths to costs, using a VSL of $10.9 million
# but we want to plot in terms of billions, so coversion is adjsuted accordingly
def deg2rad(x):
    return x * (10.9*10**-3.0)
def rad2deg(x):
    return x / (10.9*10**-3.0)
# create figure and plot datasets
fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(8,5))
ax.bar(0,tot_plot1,color='#a6cee3',yerr=[[tot_plot1-low_ci1],[hi_ci1-tot_plot1]],
       capsize=3)
ax.bar(1,tot_plot,color='#1f78b4',yerr=[[tot_plot-low_ci],[hi_ci-tot_plot]],
       capsize=3)
secax = ax.secondary_yaxis('right', functions=(deg2rad, rad2deg))
secax.set_ylabel('Savings [billions, $2019]',fontsize=16)
ax.set_xticks([0,1],['ABI-1pm','ABI-Daytime'],fontsize=16)
ax.set_ylabel('Averted Deaths\n [deaths/year]',fontsize=16)
plt.tight_layout()
plt.savefig(out_fig_path+'nat_mortalities_ervs_bar'+pct_reduce+'.png',dpi=400)
#fig.show()
plt.close()

# print numbers used in figure discussion
print('cal only totals deaths')
print(np.nansum(abi_mort[0,:,:]*cal_area_mask)-np.nansum(abi_mort_bm[0,:,:]*cal_area_mask))
print(np.nansum(abi_mort[0,:,:]*cal_area_mask)-np.nansum(abi1pm_mort_bm[0,:,:]*cal_area_mask))

#%% print additional numbers for paper
print('ERV totals')
print('baseline',np.nansum(geo_ervs[0,:,:]*area_mask))
print('geo reduced pm',np.nansum(geo_erv_bm[0,:,:]*area_mask))
print('leo reduced pm',np.nansum(leo_erv_bm[0,:,:]*area_mask))

print('cal only totals ervs')
print(np.nansum(geo_ervs[0,:,:]*cal_area_mask))
print(np.nansum(geo_erv_bm[0,:,:]*cal_area_mask))
print(np.nansum(leo_erv_bm[0,:,:]*cal_area_mask))

#%% extra plots for talks, nut used in manuscript
'''
# difference plot of code orange days in ABI1pm and ABI
diff = abi_nalerts-abi1pm_nalerts
fig, axarr = plt.subplots(nrows=1,ncols=1,subplot_kw={'projection': ccrs.PlateCarree()},
                          figsize=(8,4))
cs = axarr.pcolormesh(glon,glat,diff*plt_mask,cmap=my_magma,
                       vmin=1,vmax=25)
axarr.axis("off")
axarr.add_feature(cfeature.STATES.with_scale('50m'),edgecolor='gray')
cbar = fig.colorbar(cs,ax=axarr,orientation='horizontal',pad=0,shrink=0.6,
                    ticks = [1,5,10,15,20,25],extend='max')
cbar.set_label(label='Additional Alert Days w/ Geo',size=16)
fig.tight_layout()
plt.savefig(out_fig_path+'annual2020_alertdiffs_nat.png',dpi=300)
fig.show()
'''

# former version of fig3 with a cumulative sum
'''
us_mask = np.where(np.isnan(plt_mask),0,1)
pop_us = pop*us_mask
pop_flat = pop_us.flatten()

abi_list = np.argsort(abi_nobs.flatten())
viirs_list = np.argsort(viirs_nobs.flatten())
abi1pm_list = np.argsort(abi1pm_nobs.flatten())
an_sorted = an_df.sort_values(by='obs_count')
aqs_mon_sorted = aqs_mon.sort_values(by='obs_count')

cum_sum_plot = [np.cumsum(pop_flat.flatten()[abi_list])/np.nansum(pop_flat),
                np.cumsum(pop_flat.flatten()[viirs_list])/np.nansum(pop_flat),
                np.cumsum(pop_flat.flatten()[abi1pm_list])/np.nansum(pop_flat),
                np.nancumsum(an_sorted['pop'].values)/np.nansum(an_sorted['pop'].values),
                np.cumsum(aqs_mon_sorted['CL8AA2020'].values)/np.nansum(aqs_mon_sorted['CL8AA2020'].values)]
fig, ax = plt.subplots(1,2,figsize=(9,5))
ax[0].plot(abi_nobs.flatten()[abi_list],cum_sum_plot[0], color=colors[4])
ax[0].plot(viirs_nobs.flatten()[viirs_list],cum_sum_plot[1],color = colors[2])
ax[0].plot(abi1pm_nobs.flatten()[abi1pm_list],cum_sum_plot[2],color = colors[3])
ax[0].plot(an_sorted['obs_count'].values,cum_sum_plot[3],color = colors[1])
ax[0].plot(aqs_mon_sorted['pm_obs_count'].values,cum_sum_plot[4],color = colors[0])
# print numbers for figure discussion
inds = []
for i in range(4):
    diff = abs(cum_sum_plot[i]-0.50)
    inds.append(np.where(diff == np.min(diff))[0])
print(cum_sum_plot[0][inds[0]],abi_nobs.flatten()[abi_list[inds[0]]])
print(cum_sum_plot[1][inds[1]],viirs_nobs.flatten()[viirs_list[inds[1]]])
print(cum_sum_plot[2][inds[2]],abi1pm_nobs.flatten()[abi1pm_list[inds[2]]])
print(cum_sum_plot[3][inds[3]],an_sorted['obs_count'].values[inds[3]])
'''


