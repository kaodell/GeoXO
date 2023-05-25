#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
count_airnowalert_coverage.py
    python script to load files output from airnowequivalent code from Barron and count
    n obs and alert days by zip code
    v1 - modified to work with daily mean/max values from AirNow, output of airnow_2020_dailyareameanmax.py code from Barron
    output fn changed to '_dailymean'
Created on Wed Nov  9 09:16:58 2022
@author: kodell
"""
#%% user inputs
prj_folder = '/Users/kodell/Library/CloudStorage/GoogleDrive-kodell@email.gwu.edu/My Drive/Ongoing Projects/GeoXO/'
# zip code shapefile and population estimates from the 2020 american community survey, accessed through IPUMS NHGIS
# link: https://www.nhgis.org/
zip_shp_fn = prj_folder + 'population_data/ACS_IPUMS/2020pop_2010zctas/nhgis0011_shape/nhgis0011_shapefile_tl2020_us_zcta_2010/US_zcta_2010_tl20.shp'
zip_pop_fn = prj_folder + 'population_data/ACS_IPUMS/2020pop_2010zctas/nhgis0011_csv/nhgis0011_ts_geog2010_zcta.csv'
# zipcode to county crosswalk file to identify states in which zipcodes lie to remove duplicates from 
# airnow_2020_dailyareameanmax.py, file from 
zip_state_fn = '/Users/kodell/Desktop/uscb_zcta10_county_crosswalk.txt'#prj_folder + 'ZIP_COUNTY_032020.csv'
st_fips_fn = prj_folder + 'state_fips_usps.csv'

airnow_fp = '/Users/kodell/Desktop/airnowlike_daily/2020/ZIP_daily_'



#%% load modules
import pandas as pd
import numpy as np
import datetime as dt
import sys
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import geopandas as gpd

#%% load data and make array of dates
zip_pop = pd.read_csv(zip_pop_fn)
zip_state = pd.read_csv(zip_state_fn)
st_fips = pd.read_csv(st_fips_fn)

date = dt.date(year = 2020,month=1,day=1)
end_date = dt.date(year=2020,month=12,day=31)
date_strs = []
while date <= end_date:
    date_strs.append(date.strftime('%Y%m%d'))
    date += dt.timedelta(days=1)

#%% loop through days and count obs and alerts by zip code
# will also all population to the output file for ease of calculations later

obs = np.zeros(36713)
pm_sum = np.zeros(36713)
alerts_mean = np.zeros(36713)
alerts_max = np.zeros(36713)
for date in date_strs:
    fn_load = airnow_fp + date + '.csv'
    #airnow = pd.read_csv(fn_load)
    airnow=pd.read_csv(fn_load)
    print(airnow.shape[0])

    for i in range(airnow.shape[0]):
        if np.isnan(airnow['FINAL_mean_value'].iloc[i]):
            continue
        elif np.isnan(airnow['FINAL_max_value'].iloc[i]):
            continue
        obs[i] += 1
        pm_sum[i] += airnow['FINAL_max_value'].iloc[i] # want to use max for this
        if airnow['FINAL_mean_value'].iloc[i] >= 35.45:
            alerts_mean[i] += 1
        if airnow['FINAL_max_value'].iloc[i] >= 35.45:
            alerts_max[i] += 1

        # print(ind)
    print(date)
pm_mean = pm_sum/obs
#%% create dataframe, clean, and save

# create new dataframe
ann_airnow1 = pd.DataFrame(data = {'geoid10':airnow['GEOID10'].values,
                                  'stusps':airnow['STUSPS'].values,
                                  'obs_count':obs,
                                  'alerts_mean':alerts_mean,
                                  'alerts_max':alerts_max,
                                  'pm_mean':pm_mean})

#%% check duplicates - airnow code creates duplicates for zipcodes close to or overlapping states
# first, flag multiples
inds_multi = np.where(ann_airnow1.duplicated(subset='geoid10',keep=False))[0]
ann_airnow1['multi_flag'] = np.zeros(ann_airnow1.shape[0])
ann_airnow1.loc[inds_multi,'multi_flag'] = 1
# remove multiples with duplicated values
ann_airnow = ann_airnow1.drop_duplicates(subset = ['geoid10','pm_mean','obs_count','alerts_mean','alerts_max'],
                          keep = 'first')
ann_airnow.reset_index(inplace=True,drop=True)
dups_id = np.where(ann_airnow.duplicated(subset='geoid10',keep=False)) # keep=False flags all duplicates

# now look through duplicate zctas with different obs/alert/pm values
inds_drop = []
for dup in dups_id[0]:
    # check what state this zip code is in according to our file and the crosswalk file
    state_file = str(ann_airnow['stusps'].iloc[dup])
    state_actual_fips_list = zip_state['GEOID'].loc[zip_state['ZCTA5']==ann_airnow['geoid10'].iloc[dup]].values
    state_fips = []
    for i in range(len(state_actual_fips_list)): 
        # can be multiple matches if zipcode is in two counties
        state_fips.append(int(str(state_actual_fips_list[i]).zfill(5)[:2]))
    if len(np.unique(state_fips))==1:
        # zipcode lies in a single state, but was assinged two states in the 
        # pandas merge due to close proximity to another state, 
        # if this duplicate is not in the correct state, drop it
        state_actual_fips = np.unique(state_fips)[0]
        if state_actual_fips in [2,60,66,15,72,78]: 
            # skip these, they aren't in our st_fips file and we'll drop them later
            continue
        # get the correct state abbriviation to match state_file format
        state_actual = st_fips['usps'].loc[st_fips['fips']==state_actual_fips].values[0]
        # is this the correct state? if not, drop.
        if state_actual != state_file:
            inds_drop.append(dup)
    else:
        # here the zipcode lies in two states, and one or both is a special state, so we have different values
        # because reporting area assignment to zipcodes is treated differently in special states
        # (see header of AirNowEquivalent.pdf)
        # worth noting none of these have alert days, but one has a difference in obs of ~10 days
        print('actual multi',ann_airnow[['geoid10','alerts_max','obs_count','stusps']].iloc[dup])
        # in this case we will pick the special state
        if state_file not in ['MD','NY','MA','VT','NH','NJ']: 
            # maine is also a special state, but we need to drop one duplicate in maine and nh, so pick nh
            # this only happens for one zcta, 03579
            inds_drop.append(dup)
            
# remove duplicates
ann_airnow2 = ann_airnow.drop(inds_drop,axis=0)
ann_airnow2.reset_index(inplace=True,drop=True)

# this should be empty
dups_check = np.where(ann_airnow2.duplicated(subset='geoid10',keep=False))

#%% remove zips outside contig us
inds_drop2 = []
missing = []
for i in range(ann_airnow2.shape[0]):
    if ann_airnow2['stusps'].iloc[i] in ['AK','AS','GU','HI','PR','VI']:
        inds_drop2.append(i)
# remove these
ann_airnow3 = ann_airnow2.drop(inds_drop2,axis=0)
ann_airnow3.reset_index(inplace=True,drop=True)

#%% add population to dataframe
pop = []
for i in range(ann_airnow3.shape[0]):
    geoid = ann_airnow3['geoid10'].iloc[i]
    ind = np.where(zip_pop['ZCTAA']==geoid)[0]
    if len(ind)==0:
        print(geoid)
        pop.append(np.nan)
        continue
    pop.append(zip_pop['CL8AA2020'].iloc[ind[0]])

ann_airnow3['pop'] = pop

#%% write to csv
ann_airnow3.to_csv(prj_folder+'airnowalert_counts_dailymean.csv')


#%% plot to make sure we got em all!
# load shapefile as a geopandas array
zipdf = gpd.read_file(zip_shp_fn)
zipprojdf = zipdf.to_crs(ccrs.PlateCarree())
zipprojdf = zipprojdf.astype({"GEOID10": np.int64})
zipprojdf = zipprojdf.astype({"ZCTA5CE10": np.int64})
# merge shapefile geometries with our airnow reporting areas dataset
ann_airnow3_plot = ann_airnow3.merge(zipprojdf[['GEOID10', 'geometry','ZCTA5CE10']], 
                                     right_on='ZCTA5CE10', left_on='geoid10')
# convert airnow dataset to a geopandas array for plotting
ann_airnow3_plot1 = gpd.GeoDataFrame(ann_airnow3_plot,crs = ccrs.PlateCarree(),
                                     geometry = 'geometry')
# plot all ZCTAs in airnow dataset to make sure we got them all
fig, ax1 = plt.subplots(nrows=1,ncols=1,subplot_kw={'projection': ccrs.PlateCarree()},
                          figsize=(10,5))
ax = ann_airnow3_plot1.plot('geoid10', ax = ax1, vmin=100000,vmax=200000,cmap='magma', 
                   edgecolor='none',aspect='auto')
ax.axis("off")
ax.set_title('check for all ZCTAs/n gaps expected in remote areas and national parks',fontsize=12)
plt.show()



