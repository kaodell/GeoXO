#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
count_AQS_coverage.py
    python script to count population covered by AQS monitors using different metrics
Created on Tue Nov  1 14:41:42 2022
written by: Katelyn O'Dell, with guidance from code written by Barron H. Henderson
"""
#%% user inputs
# local path to project folder
prj_folder = '/Users/kodell/Library/CloudStorage/GoogleDrive-kodell@email.gwu.edu/My Drive/Ongoing Projects/GeoXO/'

# path to census tract shape files and population from IPUMS
ct_shp_fn = prj_folder + 'population_data/ACS_IPUMS/2020pop_2010tracts/nhgis0012_shape/nhgis0012_shapefile_tl2020_us_tract_2010/US_tract_2010_tl20.shp'
pop_data = prj_folder + 'population_data/ACS_IPUMS/2020pop_2010tracts/nhgis0012_csv/nhgis0012_ts_geog2010_tract.csv'

# url to airnow daily files
urlroot = 'https://s3-us-west-1.amazonaws.com//files.airnowtech.org/airnow'
# from this website we will load the daily_data_v2.dat files in a loop below
# documentation on these files is available here: https://docs.airnowapi.org/docs/DailyDataFactSheet.pdf

#%% import modules
import pandas as pd
import geopandas as gpd
import numpy as np
import pyproj
import sys

#%% load census tract data
gdf1 = gpd.read_file(ct_shp_fn)
pop_df = pd.read_csv(pop_data)

#%% add population to gdf
# create geoid in pop df
geoid = pop_df['STATEA'].astype('string').str.zfill(2) +\
    pop_df['COUNTYA'].astype('string').str.zfill(3)+\
    pop_df['TRACTA'].astype('string').str.zfill(6)
pop_df['geoid'] = geoid    
# merge the population and shapefile arrays
gdf = gdf1.merge(pop_df[['geoid','CL8AA2020']],left_on='GEOID10',right_on='geoid')

#%% merge mondf and add geoid
# create projection to use that matches
proj = pyproj.Proj(proj='lcc', lat_0=40, lon_0=-97, lat_1=33, lat_2=45,
                   x_0=2412000,y_0=1620000,R=6370000, 
    preserve_units=True)

# convert census tract shapefile to the same projection
gdf_use = gdf.to_crs(proj.srs)

#%% loop through days and look for observations in each census tract
# loop through dates and pull airnow data
dates = pd.date_range('2020-01-01', '2020-12-31', freq='d')

lat_missing = []
lon_missing = []

pm_obs_count_all = np.zeros(gdf.shape[0])
alert_flag_all = np.zeros(gdf.shape[0])

for date in dates:
  # create daily arrays
  pm_obs_count = np.zeros(gdf.shape[0])
  alert_flag = np.zeros(gdf.shape[0])
  # load daily file
  dvdf = pd.read_csv(f'{urlroot}/{date:%Y/%Y%m%d}/daily_data_v2.dat', delimiter='|', names='Valid date|AQSID|site name|parameter name|reporting units|value|averaging period|data source|AQI|AQICategory|latitude|longitude|fullaqsid'.split('|'))
  # pull 24hr pm2.5 data
  pm_dvdf = dvdf[dvdf['parameter name']=='PM2.5-24hr']
  pm_dvdf.reset_index(inplace=True,drop=True)
  # check for nans or very low values (shouldn't have any)
  if len(np.where(pd.isna(pm_dvdf['value']))[0]) > 0:
      sys.exit('nans in pm file')
  if len(np.where(pm_dvdf['value']<-5)[0]) > 0:
      sys.exit('possible invalid pm value')   
  # create x,y points for all monitor lat/lons in the file in the same projection as census tracts
  X, Y = proj(pm_dvdf['longitude'].values,pm_dvdf['latitude'].values)
  # use these to create a geopandas array
  geo_pm_df1 = gpd.GeoDataFrame(dict(lon=pm_dvdf['longitude'].values, lat=pm_dvdf['latitude'].values, 
                                    aqsid = pm_dvdf['AQSID'].values,value = pm_dvdf['value'].values,
                                    X=X,Y=Y),geometry=gpd.points_from_xy(X, Y), crs=proj.srs)
  # join arrays where the census tracts and points intersect
  # this adds the intersecting census tracts to the monitor array
  # documentation: https://geopandas.org/en/stable/docs/reference/api/geopandas.sjoin.html
  geo_pm_df = geo_pm_df1.sjoin(gdf_use[['GEOID10','CL8AA2020','geometry']],
                             predicate='intersects',how='left')
  # loop through new combined array to count observations and alert days for census tracts
  for i in range(geo_pm_df.shape[0]):
      geo_ind = np.where(geo_pm_df['GEOID10'].iloc[i]==gdf['GEOID10'])[0]
      if len(geo_ind)==0:
          # we're going to plot these to make sure they just arent in the US
          lat_missing.append(geo_pm_df['lat'].iloc[i])
          lon_missing.append(geo_pm_df['lon'].iloc[i]) 
      else:
          # if an observation, count the tract as covered
          pm_obs_count[geo_ind] = 1 
          # use '=' istead of '+=' because a census tract could have more than one monitor, 
          # but we want to just count the day as covered, not # of obs per day per tract
          
          # next, count alerts
          if geo_pm_df['value'].iloc[i] > 35.45:
              alert_flag[geo_ind] = 1
  pm_obs_count_all += pm_obs_count
  alert_flag_all += alert_flag
  print(f'{date:%Y-%m-%d}', end='.')
gdf['obs_count'] = pm_obs_count_all
gdf['alert_count'] = alert_flag_all

#%% write to csv
gdf.to_csv(prj_folder + 'concentration_data/EPA_AQI/AQS_obs_alert_counts_2020.csv')

#%% plots to check 
import matplotlib.pyplot as plt
import cartopy.crs as crs
import cartopy.feature as cfeature

# first, that all our missing ids are not in the US
fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(1,1,1, projection=crs.Robinson())
#ax.set_extent([-110, -125, 45, 50],crs=crs.PlateCarree())
ax.set_global()
ax.add_feature(cfeature.COASTLINE, edgecolor="tomato")
ax.add_feature(cfeature.BORDERS, edgecolor="tomato")
ax.gridlines()
plt.scatter(x=lon_missing, y=lat_missing,color="dodgerblue",s=1,alpha=0.5,transform=crs.PlateCarree()) 
plt.show()

# map of alerts in the US by census tract
gdf.plot(column='obs_count')


