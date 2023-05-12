#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
grid_baseline_mort.py
    python script to grid county-level baseline mortality rates 
    from CDC WONDER database
Created on Wed Sep 14 13:18:31 2022
written by: Katelyn O'Dell

v0 - inital code, testing
v1 - first run on Pegasus, 09.15.22
v2 - fill missing values with state mortality rates, 09.20.22 (2a = testing output)
v2_Jian - v2 code but loading 12km WRF-Chem grid from Jian He
final - final for geoXO manuscript pre-submission
"""
#%% user inputs
# overall project folder
prj_folder = '/GWSPH/home/kodell/GeoXO/datafiles/'

# CDC WONDER baseline mortality rates
# download from https://wonder.cdc.gov/
# county-level
br_fp = prj_folder + 'in/Underlying Cause of Death, 1999-2020.txt'
# state-level for filling in missing counties
state_br_fp = prj_folder + 'in/Underlying Cause of Death, 1999-2020 State.txt'

# grid for satellite-based pm2.5 datasets to regrid the rates to
to_grid_fp = '/GWSPH/home/kodell/GeoXO/datafiles/in/ABI_daily/2020/ABI_CONUS_PM2.5_daily_20200101_0.01x0.01.nc'

# 2019 TIGER/Line shapefiles for counties
# downlaod from https://www2.census.gov/geo/tiger/TIGER2019/COUNTY/
shp_fp = prj_folder + 'in/tl_2019_us_county/tl_2019_us_county.shp'

# indicate where to put output data
out_fp = prj_folder + 'out/'

# verison to label out files, see header for version descriptions
version = 'final'

#%% load modules
import pandas as pd
import geopandas as gpd
from netCDF4 import Dataset
import numpy as np

#%% import data
# county shape file
cnty_shp_df1 = gpd.read_file(shp_fp)

# grid to put gridded mortality rates on
nc_fid = Dataset(to_grid_fp)
lat = nc_fid['lat'][:]
lon = nc_fid['lon'][:]
nc_fid.close()
glon1, glat1 = np.meshgrid(lon,lat)
glat = glat1.flatten()
glon = glon1.flatten()

# baseline rates files
# county level, load rows before totals and file information at the end of the csv
br_df = pd.read_csv(br_fp,sep='\t',nrows=3147,dtype='str')
# state level, load rows before totals and file information at the end of the csv
state_br_df = pd.read_csv(state_br_fp,sep='\t',nrows=51,dtype='str')

#%% crosswalk county shapefile and baseline rates file

# update FIPS 46113 (Shannon cnty SD) which was renamed to Oglala Lakota County (FIPS 46102)
ind_up = np.where(br_df['County Code'] == '46113')[0]
br_df.loc[ind_up,'County Code'] = '46102'

# drop values outside the contig US
codes_drop = ['02','15','60','66','69','72','78']
for code_drop in codes_drop:
    cnty_shp_df1.drop(cnty_shp_df1[cnty_shp_df1['STATEFP']==code_drop].index,inplace=True)
cnty_shp_df1.reset_index(inplace=True,drop=True)

# first combine the county shapefile dataframe with the baseline rate dataframe
# using pandas merge function (documentation: https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.merge.html)
comb = cnty_shp_df1.merge(br_df,left_on='GEOID',right_on='County Code')

# rename crude rate the in the state file, then merge with the combined dataframe
state_br_df['Crude Rate State'] = state_br_df['Crude Rate'].copy(deep='True')
comb2 = comb.merge(state_br_df[['State Code','Crude Rate State']],
                   left_on='STATEFP',right_on='State Code')
# find inds where county-level mort is labeled as unreliable or suppressed
unreliable_inds = np.where(comb2['Crude Rate']=='Unreliable')[0]
suppressed_inds = np.where(comb2['Crude Rate']=='Suppressed')[0]
# create a final version of crude rates where these missing values are replaced with state level mort rate
comb2['Crude Rate Final'] = comb2['Crude Rate'].copy(deep='True')
comb2.loc[unreliable_inds,'Crude Rate Final']= comb2['Crude Rate State'].iloc[unreliable_inds].values
comb2.loc[suppressed_inds,'Crude Rate Final'] = comb2['Crude Rate State'].iloc[suppressed_inds].values
comb2['Crude Rate Final']=comb2['Crude Rate Final'].astype('float')
cnty_shp_df = comb2.copy(deep=True)

# find counties in baseline rates but not in shapefile
# all alaska or HI, which we removed from the shapefile
# and FIPS 46113 (Shannon cnty SD) which was renamed to Oglala Lakota County (FIPS 46102)
# there are also two independent cities in virginia (Bedford and Clifton Forge)
# that later merged into counties, but their rates are labelled as "missing" so we can't use them anyway
# fixed the FIPS change above by updating the fips code in the baseline rates
'''
for bi in range(br_df.shape[0]):
    cnty_id = str(br_df['County Code'].iloc[bi]).zfill(5)
    ind = np.where(cnty_shp_df['GEOID']==cnty_id)[0]
    if len(ind)<1:
        print(cnty_id)
'''

#%% create shaply points from lat/lon coordinates using geopandas
points_df = pd.DataFrame(data = {'lat':glat,'lon':glon})
points_gdf1 = gpd.GeoDataFrame(points_df,geometry=gpd.points_from_xy(points_df.lon,points_df.lat))
points_gdf = points_gdf1.set_crs(cnty_shp_df1.crs) # match projection of TIGER/line shapefiles

#%% loop through county and assign county baseline value to grid cell
gmort = np.empty(glon.shape)
gmort[:] = -7777.7
# also test that the grid cells are correctly assinged
gtest_mort = np.zeros(glon.shape)
for i in range(cnty_shp_df.shape[0]):
    # find which points are inside county i
    ind = np.where(points_gdf.within(cnty_shp_df['geometry'].iloc[i]).values)
    # assign these values the county mortality rate
    gmort[ind]=cnty_shp_df['Crude Rate Final'].iloc[i]
    # count which grid cells are assigned
    gtest_mort[ind] += 1 # check to see if any grid cells are double assigned, and can see which grid cells aren't assinged
    print(i)
# add this data to the dataframe
points_df['cnty_mort'] = gmort

#%% save data
gmort_df = pd.DataFrame(data = {'lat':glat,'lon':glon,'cnty_mort_rate':gmort,'assignment_counts':gtest_mort})
gmort_df.to_csv(out_fp + 'gridded_cnty_mort_'+version+'.csv')

# we will plot this and check in the plot_national_results.py code, which is run locally
