#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
calc_annual_alert_HIA.py
    v1. python script to read and prepare shobha's data to use in our analysis
    v2. updated to work with the higher resolution data, 0.01
        output names changed to '0.01' to match data inputs & no longer need to 
        combine the goes east/west grid
    v3. editied code to run mortality HIA and annual averages for viirs and abi 
    v3b. testing for version 4
    v4. fixed popw and areaw averaging and changed mort baseline to national
    separately, changed name to calc_annual_alert_HIA.py
    v5. reorganized structure to test sensitivity cases: proxy, og, sat, and abi-only
    v6. switch to using 1pm ABI for the VIIRS case
    v7. include both the VIIRS and ABI cases, reorganize code
    v8. add county-level baseline mortality rates
    v9. updated 07/22 ABI, daily and ABI, 1pm files
    v10. final pre-submission
Created on Tue Feb  1 13:44:19 2022
@author: kodell
"""
#%% user inputs
# run options
# three arguments passed in the submission script:
    # 1) case identifies which dataset to use: abi, abi1pm, or viirs
    # 2) case_desc indicates to run original or proxy (substitute abi values)
    #   abi1pm should use proxy and viirs should use og to match the manuscript
    # 3) pct_reduce_str is a string indicating the percent reduction in PM for behavior modification.
    #   0.3 is used in the main text; 0.15 and 0.5 are run for supplemental sensitivity analysis
import sys, os
case = sys.argv[1]
case_desc = sys.argv[2]
pct_reduce_str = sys.argv[3]

# file locations
# folder which contains pm2.5 datasets
prj_folder = '/GWSPH/home/kodell/GeoXO/'
# specific dataset paths
viirs_pm_file = prj_folder + 'datafiles/in/VIIRS/2020/VIIRS_CONUS_PM2.5_daily_' 
abi1pm_pm_file = prj_folder + 'datafiles/in/ABI1pm/2020/ABI1pm_CONUS_PM2.5_daily_' 
abi_pm_file = prj_folder + 'datafiles/in/ABI_daily/2020/ABI_CONUS_PM2.5_daily_' 
abi1pm_files = os.listdir(prj_folder + 'datafiles/in/ABI1pm/2020/')
# regridded population file
pop_file = prj_folder + 'datafiles/in/population_data/regrided_2020pop_0.01_final.nc'
# gridded baseline mortality rates
br_file = prj_folder + 'datafiles/out/gridded_mort/gridded_cnty_mort_final.csv'

# description and path of output file 
out_desc = '0.01_'+case+'_'+case_desc+'_v10_'+pct_reduce_str 
out_fp = prj_folder + 'datafiles/out/'

# alert day cutoffs and percent reduction in pm exposure
alert_cutoff = 35.45
alert_cutoff_aqi = 101.0
pct_reduce = float(pct_reduce_str)

# hard code relative risks for all-cause mortality from long-term exposure
mort_RRs = [1.06,1.04,1.08] # Turner et al., 2016; Table E10; single pollutant model, same as Castillo et al., 2021

#%% import modeules
from netCDF4 import Dataset
import numpy as np
import datetime as dt
import pandas as pd

#%% user-defined functions
def chronic_HIA(conc, cf, pop, base_rate, betas, grid_area): 
    # conc = concentration
    # cf = counter-factual concentration, or level at which no excess risk is assumed
    # pop = population on the same grid as concentration
    # base_rate = baseline mortality rate on the same grid as concentration 
    # betas = array of beta values: [beta_central, beta_lower, beta_upper] where lower and upper
    #           are the lower and upper bounds of the 95% confidence iterval, respectively
    # grid_area = area of grid cells for the concentration grid
    paf_avg_out = [] # population attributable fraction in percent
    events_tot_out = [] # attributable mortality
    events_tot_pp_out = [] # attributable mortality per 10^5 people
    events_tot_pk_out = [] # attributable mortality per kilometer
    
    z = np.where(conc<cf,0,conc-cf)
    for beta in betas:      
        paf_avg = 100.0*(1.0 - np.exp(-beta*z))
        events = (paf_avg/100.0)*pop*(base_rate)
        events_tot = events
        events_tot_pk = (1000/grid_area)*events_tot
        events_tot_pp = events_tot/(pop/100000)
        
        paf_avg_out.append(paf_avg)
        events_tot_out.append(events_tot)
        events_tot_pp_out.append(events_tot_pp)
        events_tot_pk_out.append(events_tot_pk)          
    return paf_avg_out, events_tot_out, events_tot_pp_out, events_tot_pk_out

#%% load data
# gridded population counts
nc_fid = Dataset(pop_file)
pop = nc_fid['population'][:].data
grid_area = nc_fid['grid_area'][:].data
area_mask = nc_fid['us_mask'][:].data
lat = nc_fid['lat'][:].data
lon = nc_fid['lon'][:].data
nc_fid.close()

# county-level baseline mortality rate
gbr_df = pd.read_csv(br_file)
br_mort = gbr_df['cnty_mort_rate'].values
br_mort_grid1 = br_mort.reshape([2600, 5900])
br_mort_grid = br_mort_grid1/100000.0

#%% make array of dates to load files
start_day = dt.date(year=2020,month=1,day=1)
end_day = dt.date(year=2020,month=12,day=31)
dates = [start_day]
date_str = [start_day.strftime('%Y%m%d')]
day = start_day
while day < end_day:
    day = dates[-1]+dt.timedelta(days=1)
    dates.append(day)
    date_str.append(day.strftime('%Y%m%d'))
    
#%% load data and create reduced pm for selected case
# base satellite-based pm25
sat_pm25 = np.empty([len(dates),lat.shape[0],lon.shape[0]])
# behavior modification reduced pm25 using epa percent reduction values
pm25_bm24_epa = np.empty([len(dates),lat.shape[0],lon.shape[0]])
di = 0
for date in date_str:    
    abi = Dataset(abi_pm_file + date +'_0.01x0.01.nc')
    abi_day = abi['pm25'][:].data    
    abi_day = np.where(abi_day==-9999.0,np.nan,abi_day)
    abi.close()
    if case in ['abi1pm','viirs']:
        if case == 'abi1pm': # check for missing days - 0603 is missing 1pm obs per email from Zigang on 08/17/22
            pm_fp = abi1pm_pm_file + date +'_0.01x0.01.nc'
            pm_fn = pm_fp.split('/')[-1]
            if pm_fn in abi1pm_files:
                pm = Dataset(pm_fp)
                pm_day = pm['pm25'][:].data
                pm.close()
            else:
                pm_day = np.empty([pop.shape[0],pop.shape[1]])
                pm_day[:] = np.nan
                print(date,'missing')

        elif case == 'viirs':
            pm_fp = viirs_pm_file + date +'_0.01x0.01.nc'
            pm = Dataset(pm_fp)
            pm_day = pm['pm25'][:].data
            pm.close()
        pm_day = np.where(pm_day==-9999.0,np.nan,pm_day)

        # now loop through the two case descriptions
        # proxy = subsititue ABI values for the array
        # og = use original values for the array
        if case_desc == 'proxy':
            pm_day_use = np.where(np.isnan(pm_day),np.nan,abi_day)
            sat_pm25[di,:,:] = pm_day_use
            pm25_bm24_epa[di,:,:] = np.where(pm_day_use>=alert_cutoff, pm_day_use*(1-pct_reduce),abi_day)
        elif case_desc == 'og':
            sat_pm25[di,:,:] = pm_day
            pm25_bm24_epa[di,:,:] = np.where(pm_day>=alert_cutoff, pm_day*(1-pct_reduce),pm_day)
        else:
            sys.exit('case description not understood') 
        del abi_day, pm_day
        
    elif case =='abi':
        sat_pm25[di,:,:] = abi_day
        pm25_bm24_epa[di,:,:] = np.where(abi_day>=alert_cutoff, abi_day*(1-pct_reduce),abi_day)
        del abi_day

    else:
        sys.exit('case not found')

    di += 1
            
#%% calcualte annual averages and number of observations
ann_avg_PM = np.nanmean(sat_pm25, axis=0)
ann_avg_pm_bm24_epa = np.nanmean(pm25_bm24_epa,axis=0)
ann_count = np.sum(np.where(np.isnan(sat_pm25),0,1),axis=0)

#%% run mortality  HIA
mort_betas = np.log(mort_RRs)/10.0

# calculate reduction in mortality from this reduction in long-term exposure
[paf_avg_out, base_mort_tot_out, events_tot_pp_out, events_tot_pk_out] = chronic_HIA(ann_avg_PM, 2.8, pop, 
                                                                                     br_mort_grid,mort_betas, grid_area)
[paf_avg_out, BM_mort_tot_out, events_tot_pp_out, events_tot_pk_out] = chronic_HIA(ann_avg_pm_bm24_epa, 2.8, pop, 
                                                                                      br_mort_grid, mort_betas, grid_area)
# reshape these arrays
base_mort_tot_out2 = np.reshape(base_mort_tot_out,[3,base_mort_tot_out[0].shape[0],base_mort_tot_out[0].shape[1]])
BM_mort_tot_out2 = np.reshape(BM_mort_tot_out,[3,BM_mort_tot_out[0].shape[0],BM_mort_tot_out[0].shape[1]])

#%% calculate additional annual averages/data

# 2) population-weighted daily PM
pm_mask = np.where(np.isnan(sat_pm25),0,1)
pm_mask_bm = np.where(np.isnan(pm25_bm24_epa),0,1)
PM_popw = np.nansum(sat_pm25*pop*area_mask,axis=(1,2))/np.nansum(pop*area_mask*pm_mask,axis=(1,2))
PM_bm_popw = np.nansum(pm25_bm24_epa*pop*area_mask,axis=(1,2))/np.nansum(pop*area_mask*pm_mask_bm,axis=(1,2))

# 3) area-weighted daily PM
PM_areaw = np.nansum(sat_pm25*grid_area*area_mask,axis=(1,2))/np.nansum(grid_area*area_mask*pm_mask,axis=(1,2))
PM_bm_areaw = np.nansum(pm25_bm24_epa*grid_area*area_mask,axis=(1,2))/np.nansum(grid_area*area_mask*pm_mask_bm,axis=(1,2))

# 4) total coverage and alert days for each dataset
PM_obs = np.sum((pm_mask*area_mask),axis=0)
PM_alerts = np.sum(np.where((sat_pm25*area_mask)>=alert_cutoff,1,0),axis=0)

# 5) mean and standard deviation of spatial coverage
daily_area_covered = np.nansum(pm_mask*grid_area*area_mask,axis=(1,2))

#%% save out combined data and annual averages to use in future codes
nc_w_fid = Dataset(out_fp+'alerts_chronic_HIA_'+out_desc+'.nc', 
                           'w', clobber=True,  format='NETCDF4')
nc_w_fid.description = 'alerts HIA w/ alert cutoff:'+str(alert_cutoff)+' and '+str(pct_reduce)+'% PM reduction'
nc_w_fid.history = 'Created' + dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

nc_w_fid.createDimension('time', PM_popw.shape[0])
nc_w_fid.createDimension('gridx', lat.shape[0])
nc_w_fid.createDimension('gridy', lon.shape[0])
nc_w_fid.createDimension('uncertainty',None)
# lats and lons
lat_w = nc_w_fid.createVariable('lat', np.float32, ('gridx',))
lon_w = nc_w_fid.createVariable('lon', np.float32, ('gridy',))
# daily PM, pop weighted
daily_PM_popw_w = nc_w_fid.createVariable('daily_PM_popw', np.float32, ('time',))
daily_PMbm_popw_w = nc_w_fid.createVariable('daily_PMbm_popw', np.float32, ('time',))
# daily PM, area weighted
daily_PM_areaw_w = nc_w_fid.createVariable('daily_PM_areaw', np.float32, ('time',))
daily_PMbm_areaw_w = nc_w_fid.createVariable('daily_PMbm_areaw', np.float32, ('time',))
# coverage and alert days
PM_obs_w = nc_w_fid.createVariable('PM_obs', np.float32, ('gridx','gridy',))
PM_alerts_w = nc_w_fid.createVariable('PM_alerts', np.float32, ('gridx','gridy',))
# HIA totals - mortalities
base_mort_w = nc_w_fid.createVariable('base_mort', np.float32, ('uncertainty','gridx','gridy',))
bm_mort_w = nc_w_fid.createVariable('bm_mort', np.float32, ('uncertainty','gridx','gridy',))
# annual average pm
annavg_PM_w = nc_w_fid.createVariable('annavg_PM', np.float32, ('gridx','gridy',))
annavg_PM_bm_w = nc_w_fid.createVariable('annavg_PM_bm', np.float32, ('gridx','gridy',))
# pm mask for each day
area_covered_w = nc_w_fid.createVariable('area_covered',int,('time'))

# now put the data in the variables
# lats and lons
lon_w[:] = lon
lat_w[:] = lat
# daily PM, pop weighted
daily_PM_popw_w[:] = PM_popw 
daily_PMbm_popw_w[:] = PM_bm_popw
# daily PM, area weighted
daily_PM_areaw_w[:] = PM_areaw 
daily_PMbm_areaw_w[:] = PM_bm_areaw 
# coverage and alert days
PM_obs_w[:] = PM_obs
PM_alerts_w[:] = PM_alerts
# HIA totals - mortalities
base_mort_w[:] = base_mort_tot_out2 
bm_mort_w[:] = BM_mort_tot_out2
# annual averages
annavg_PM_w[:] = ann_avg_PM
annavg_PM_bm_w[:] = ann_avg_pm_bm24_epa
# pm area covered
area_covered_w[:] = daily_area_covered
nc_w_fid.close()








