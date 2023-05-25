#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
calc_national_ervs.py
    python script to calculate averted asthma emergency room visits with geo v leo forecasting
    for of partciulate pollution events nationally
    
    updated to run on server and save annual totals and national timeseries
    in a netCDF file
    
    updated to load daily files individually and process (v2)
    
    edited to use 1pm ABI data as  leo
    
    revised to retain 1pm ABI and leo
    
    v3. is the higher resolution data
    v4. includes airnow 
    v6. updating for ABI 07/22 processing issue from Shobha's team 
    v7. final version for PNAS submission, code title changed to calc_national_ervs.py
    
Created on Tue Mar  1 15:05:52 2022
written by: Katelyn O'Dell
"""
#%% user inputs
# file locations
# folder which contains pm2.5 datasets
prj_folder = '/GWSPH/home/kodell/GeoXO/'
# specific dataset paths
abi1pm_file = prj_folder + 'datafiles/in/ABI1pm/2020/ABI1pm_CONUS_PM2.5_daily_' 
abi_pm_file = prj_folder + 'datafiles/in/ABI_daily/2020/ABI_CONUS_PM2.5_daily_' 
import os
abi1_files = os.listdir(prj_folder + 'datafiles/in/ABI1pm/2020/')
# regridded population file
pop_file = prj_folder + 'datafiles/in/population_data/regrided_2020pop_0.01_final.nc'

# hard-code asthma baseline rates for short-term exposure analysis
# these are from HCUP via the CDC webpage on asthma in the US, values for 2018.
# https://www.cdc.gov/asthma/healthcare-use/2018/table_a.html
# last accesed 11 May 2023
asth_erv_br = 50.2/10000 # SE 2.34
# relative risk values from Orellano et al. 2017 meta analysis doi:https://doi.org/10.1371/journal.pone.0174050
asth_rrs = [1.03, 1.01, 1.05]

# set threshold for alerts
alert_cutoff = 35.45
alert_cutoff_aqi = 101.0
pct_reduce = 0.30

# case for flagging alert days
case_desc = 'all' # use proxy data for flagging

# where to put output files
version = 'v7'
out_fp = prj_folder + 'datafiles/out/'
out_desc = str(pct_reduce)[2:] + '_' + case_desc +'_'+ version

#%% import modules
from netCDF4 import Dataset
import numpy as np
import datetime as dt

#%% define acute HIA function
def acute_HIA(conc, cf, pop, base_rate, betas, area_mask):
    """Function calculates a health impact assessment for daily-average exposure,
    assuming a log normal concentration response function.
    
    Parameters
    ----------
    conc: numpy.ndarray
        array of daily-average pollutant concentrations in same unit as beta
    cf: float
        counter-factual concentration below which no health impact is assumed
    pop: numpy.ndarray
        population in each grid cell. on the same grid as conc input.
    base_rate: float
        annual baseline rate for health outcome
    betas: tuple
        floats of betas to use in the HIA.
    area_mask: numpy.ndarray
        array of 1s and 0s indicating the grid cells over which to calcaute the HIA.
        1 = include, 0 = do not include
    
    Returns
    -------
    paf_avg_out: tuple of numpy.ndarrays for each beta input
        average population attributable fraction over the time period
    events_tot_out: tuple of numpy.ndarrays for each beta input
        total attributable health events each day
    events_tot_pp_out: tuple of numpy.ndarrays for each beta input
        daily events per person in each grid cell
    events_tot_pk_out: tuple of numpy.ndarrays for each beta input
        daily events per km2 in each grid cell
    
    ------
    written by Katelyn O'Dell
    """

    events_tot_out = []
    z = np.where(conc<cf,0,conc-cf)
    for beta in betas:
        paf = 100.0*(1.0 - np.exp(-beta*z))
        events = (paf/100.0)*pop*((base_rate/365)) # calculate at a daily level
        events_tot = events
        events_tot_out.append(events_tot)
            
    return events_tot_out*area_mask

#%% load population data
# gridded population counts
nc_fid = Dataset(pop_file)
pop = nc_fid['population'][:].data
grid_area = nc_fid['grid_area'][:].data
area_mask = nc_fid['us_mask'][:].data
lat = nc_fid['lat'][:].data
lon = nc_fid['lon'][:].data
nc_fid.close()

#%% load data 
start_day = dt.date(year=2020,month=1,day=1)
end_day = dt.date(year=2020,month=12,day=31)
dates = [start_day]
date_str = [start_day.strftime('%Y%m%d')]
day = start_day
while day < end_day:
    day = dates[-1]+dt.timedelta(days=1)
    dates.append(day)
    date_str.append(day.strftime('%Y%m%d'))

# pre-allocate arrays
# total baseline ervs
geo_events_sum = np.zeros([3,pop.shape[0],pop.shape[1]])
abi1_events_sum = np.zeros([3,pop.shape[0],pop.shape[1]])
# ervs under reduced pm exposure with behavior modification
geo_events_bm_sum = np.zeros([3,pop.shape[0],pop.shape[1]])
abi1_events_bm_sum = np.zeros([3,pop.shape[0],pop.shape[1]])
# loop through dates and calculate averted ervs
di = 0
for date in date_str:
    # load abi 1pm
    abi1_fp = abi1pm_file + date +'_0.01x0.01.nc'
    abi1_fn = abi1_fp.split('/')[-1]
    # check for file first, abi 1pm is missing a day
    if abi1_fn in abi1_files:
        abi1 = Dataset(abi1pm_file + date +'_0.01x0.01.nc')
        abi1_pm25 = abi1['pm25'][:].data
        abi1.close()
    else:
        abi1_pm25 = np.empty([pop.shape[0],pop.shape[1]])
        abi1_pm25[:] = np.nan
        print(date,'missing')

    # load abi
    abi = Dataset(abi_pm_file + date +'_0.01x0.01.nc')
    abi_pm25 = abi['pm25'][:].data
    abi.close()
    
    # replace fill values with nans, and create abi1pm proxy
    abi_pm25_1 = np.where(abi_pm25==-9999.0,np.nan,abi_pm25)
    abi1_pm25_1 = np.where(abi1_pm25==-9999.0,np.nan,abi_pm25_1) # here we want to match abi values
        
    # mask values outside the US
    abi1_pm_us = abi1_pm25_1*area_mask
    abi_pm_us = abi_pm25_1*area_mask
        
    # calcuate reduced pm25 where there are alerts
    geo_pm25_bm24_epa = np.where(abi_pm_us>=alert_cutoff, (1.0-pct_reduce)*abi_pm_us ,abi_pm_us)
    abi1_pm25_bm24_epa = np.where(abi1_pm_us>=alert_cutoff, (1.0-pct_reduce)*abi1_pm_us, abi_pm_us) # need to put abi here not abi1 so the baseline matches

    # calculate betas from RRs for the HIA equation
    asth_rrs = np.array(asth_rrs)
    betas = np.log(asth_rrs)/10.0
        
    ### ABI ####
    # no alerts
    geo_events_tot =  acute_HIA(abi_pm_us, 0, pop, asth_erv_br, betas, area_mask)
    # EPA-based exposure reduction with alerts
    geo_events_tot_bm24_epa =  acute_HIA(geo_pm25_bm24_epa, 0, pop, asth_erv_br, betas, area_mask)
    
    ### ABI, 1pm ####
    # no alerts
    abi1_events_tot =  acute_HIA(abi1_pm_us, 0, pop, asth_erv_br, betas, area_mask)
    # EPA-based exposure reduction with alerts
    abi1_events_tot_bm24_epa =  acute_HIA(abi1_pm25_bm24_epa, 0, pop, asth_erv_br, betas, area_mask)
        
    # reshape these lists into numpy arrays
    geo_events_tot2 = np.reshape(geo_events_tot,[3,geo_events_tot[0].shape[0],
                                                 geo_events_tot[0].shape[1]])
    geo_events_tot_bm24_epa2 = np.reshape(geo_events_tot_bm24_epa,
                                          [3,geo_events_tot_bm24_epa[0].shape[0],
                                           geo_events_tot_bm24_epa[0].shape[1]])

    abi1_events_tot2 = np.reshape(abi1_events_tot,[3,abi1_events_tot[0].shape[0],
                                                 abi1_events_tot[0].shape[1]])
    abi1_events_tot_bm24_epa2 = np.reshape(abi1_events_tot_bm24_epa,
                                          [3,abi1_events_tot_bm24_epa[0].shape[0],
                                           abi1_events_tot_bm24_epa[0].shape[1]])

    # sum asth ED totals for the year, making nans zero
    geo_events_sum += np.where(np.isnan(geo_events_tot2),0,geo_events_tot2)
    geo_events_bm_sum += np.where(np.isnan(geo_events_tot_bm24_epa2),0,geo_events_tot_bm24_epa2)

    abi1_events_sum += np.where(np.isnan(abi1_events_tot2),0,abi1_events_tot2)
    abi1_events_bm_sum += np.where(np.isnan(abi1_events_tot_bm24_epa2),0,abi1_events_tot_bm24_epa2)
     
    # now on to next day
    di += 1
    print(date)

#%% write to netCDF
nc_w_fid = Dataset(out_fp+'national_alerts_HIA_'+out_desc+'.nc', 
                           'w', clobber=True,  format='NETCDF4')
nc_w_fid.description = 'alerts HIA w/ alert cutoff:'+str(alert_cutoff)+' and '+str(pct_reduce)+'% PM reduction'
nc_w_fid.history = 'Created' + dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

nc_w_fid.createDimension('time', len(dates))
nc_w_fid.createDimension('gridx', lat.shape[0])
nc_w_fid.createDimension('gridy', lon.shape[0])
nc_w_fid.createDimension('uncertainty',None)
nc_w_fid.createDimension('counts',None)

# lats and lons
lat_w = nc_w_fid.createVariable('lat', np.float32, ('gridx',))
lon_w = nc_w_fid.createVariable('lon', np.float32, ('gridy',))
# HIA totals - emergency room visits
geo_erv_w = nc_w_fid.createVariable('geo_erv', np.float32, ('uncertainty','gridx','gridy'))
geo_erv_bm_w = nc_w_fid.createVariable('geo_erv_bm', np.float32, ('uncertainty','gridx','gridy'))

abi1_erv_w = nc_w_fid.createVariable('abi1_erv', np.float32, ('uncertainty','gridx','gridy'))
abi1_erv_bm_w = nc_w_fid.createVariable('abi1_erv_bm', np.float32, ('uncertainty','gridx','gridy'))

obs_abi_w = nc_w_fid.createVariable('obs_abi', np.float32, ('counts'))
obs_both_w = nc_w_fid.createVariable('obs_both', np.float32, ('counts'))
obs_abi1pm_w = nc_w_fid.createVariable('obs_abi1pm', np.float32, ('counts'))

# now put the data in the variables
# lats and lons
lon_w[:] = lon
lat_w[:] = lat
# HIA totals - emergency room visits
geo_erv_w[:] = geo_events_sum
geo_erv_bm_w[:] = geo_events_bm_sum
abi1_erv_w[:] = abi1_events_sum
abi1_erv_bm_w[:] = abi1_events_bm_sum
nc_w_fid.close()



