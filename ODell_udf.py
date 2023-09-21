#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ODell_udf.py
    python script of functions written by me or by others passed on to me
Created on Wed Sep  8 09:09:22 2021
@author: kodell
"""
#%% packages needed
import numpy as np
import cartopy.feature as cfeature
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import matplotlib.pyplot as plt
import matplotlib as mplt
from matplotlib import colors
mplt.rcParams['font.size'] = '14'
mplt.rcParams['font.family'] = 'sans-serif'
#mplt.rcParams['font.sans-serif'] = 'Veranda'

#%% make a basic map of the US using cartopy
# NOTE this assumes a PlateCaree projection
def plt_map(dlon,dlat,data,size,cmap,clabel,title,**kwargs):
    """Function creates a scatter point figure using cartopy mapping modules. 
    Can create a multi-panel figure and write figure to a file given user keyword args.
    
    Parameters
    ----------
    dlon : numpy.ndarray 
        2D array of longitude in degrees
    dlat : numpy.ndarray 
        2D array of latitude in degrees
    data : numpy.ndarry or tuple of numpy.ndarrys if multi is indicated
        2D array of the data to be plotted. must be on the same grid as other
        data if indicating multiple data grids to plot
    size: size for scatter points
    cmap: tuple
        string indicating cmap to use
    clabel: tuple
        string for colorbar label
    title: tuple
        string for figure title
    
    Keyword Arguments
    -------
    clim: tuple
        interger or float limits for the colorbar colors
        otherwise limits are chosen by matplotlib defaults
    outname: string
        filepath and name for file to save figure to
        otherwise figure is not written to a file
    cpts: tuple
        interter or float of min, middle, and max of colorbar
        otherwise limits are chosen by matplotlib defaults
    multi: tuple
        [n rows, n cols] for a multi panel figure. must match shape of data,
        cmap, clabel, and title. default is [1,1].
    bkcolor: string 
        color for plot background. default is white.
    norm: string
        normalization for colorbar if desired. cannot be used with cpts.
    
    Returns
    -------
    none. figure is shown and saved if outname keyword is used.
    ------
    written by Katelyn O'Dell
    """
    vlim = kwargs.get('clim', None)
    outpath = kwargs.get('outname',None)
    vpts = kwargs.get('cpts',None)
    multi = kwargs.get('multi',None)
    bkcolor = kwargs.get('bkcolor',None)
    norm = kwargs.get('norm',None)
    shp = kwargs.get('shpfile',None)
    if multi:
        nd = len(data)
        if bkcolor:
            fig, axarr = plt.subplots(nrows=multi[0],ncols=multi[1],subplot_kw={'projection': ccrs.PlateCarree()},
                                      figsize=(11,8.5),facecolor=bkcolor)
        else:
            fig, axarr = plt.subplots(nrows=multi[0],ncols=multi[1],subplot_kw={'projection': ccrs.PlateCarree()},
                                      figsize=(11,8.5))

        axarr = axarr.flatten()
        for di in range(nd):
            ax = axarr[di]
            ax.patch.set_visible(False)
            # plot shapfile with colors
            #ax.add_feature(cfeature.LAND.with_scale('50m'),facecolor='gray',alpha=0.5)
            #ax.add_feature(cfeature.OCEAN.with_scale('50m'))

            if shp:
                reader = shpreader.Reader(shp)
                counties = list(reader.geometries())
                COUNTIES = cfeature.ShapelyFeature(counties, ccrs.PlateCarree())
                ax.add_feature(COUNTIES, facecolor='none', edgecolor='black',alpha=0.5)
            else:
                ax.add_feature(cfeature.STATES.with_scale('50m'),edgecolor='gray')
                ax.set_extent([-125, -66, 25, 50], crs=ccrs.PlateCarree())


            if bkcolor:
                ax.outline_patch.set_edgecolor(bkcolor)
            else:
                ax.outline_patch.set_edgecolor('white')
            if vlim:
                cs = ax.pcolormesh(dlon,dlat,data[di],shading='nearest',
                            transform=ccrs.PlateCarree(),cmap=cmap[di],vmin=vlim[di][0],vmax=vlim[di][1])
            elif vpts:
                divnorm=colors.TwoSlopeNorm(vmin=vpts[di][0], vcenter=vpts[di][1], vmax=vpts[di][2])
                cs = ax.pcolormesh(dlon,dlat,data[di],shading='nearest',
                            transform=ccrs.PlateCarree(),cmap=cmap[di],norm=divnorm)
            elif norm:
                cs = ax.pcolormesh(dlon,dlat,data[di],shading='nearest',
                            transform=ccrs.PlateCarree(),cmap=cmap[di],norm=norm)
            else:
                cs = ax.pcolormesh(dlon,dlat,data[di],shading='nearest',
                            transform=ccrs.PlateCarree(),cmap=cmap[di])
            cbar = fig.colorbar(cs,ax=ax,orientation='horizontal',pad=0,shrink=0.6)
            #cbar = fig.colorbar(cs,ax=ax,orientation='vertical',pad=0,shrink=0.5)
            if bkcolor:
                cbar.set_label(label=clabel[di],size=16,color='white')
                cbar.ax.xaxis.set_tick_params(color='white', labelcolor='white')
                ax.set_title(title[di],fontsize=18,color='white')
            else:
                cbar.set_label(label=clabel[di],size=16)
                ax.set_title(title[di],fontsize=18)
            plt.tight_layout()
    else:
        if bkcolor:
            fig, ax = plt.subplots(nrows=1,ncols=1,
                                      subplot_kw={'projection': ccrs.PlateCarree()},
                                      figsize=(11,8.5),facecolor=bkcolor)
            ax.outline_patch.set_edgecolor(bkcolor)

        else:
            fig, ax = plt.subplots(nrows=1,ncols=1,
                                      subplot_kw={'projection': ccrs.PlateCarree()},
                                      figsize=(11,8.5))
            ax.outline_patch.set_edgecolor('white')

        ax.patch.set_visible(False)
        # plot shapfile with colors
        #ax.add_feature(cfeature.LAND.with_scale('50m'),facecolor='white',alpha=0.5)
        #ax.add_feature(cfeature.OCEAN.with_scale('50m'))
        if shp:
            reader = shpreader.Reader(shp)
            counties = list(reader.geometries())
            COUNTIES = cfeature.ShapelyFeature(counties, ccrs.PlateCarree())
            ax.add_feature(COUNTIES, facecolor='none', edgecolor='black',alpha=0.5)
        else:
            ax.add_feature(cfeature.STATES.with_scale('50m'),edgecolor='gray')
            ax.outline_patch.set_edgecolor('white')
            ax.set_extent([-125, -66, 24, 50], crs=ccrs.PlateCarree())

        if vlim:
            cs = ax.pcolormesh(dlon,dlat,data,shading='nearest',
                        transform=ccrs.PlateCarree(),cmap=cmap,vmin=vlim[0],vmax=vlim[1])
        elif vpts:
            divnorm=colors.TwoSlopeNorm(vmin=vpts[0], vcenter=vpts[1], vmax=vpts[2])
            cs = ax.pcolormesh(dlon,dlat,data,shading='nearest',
                        transform=ccrs.PlateCarree(),cmap=cmap,norm=divnorm)
        elif norm:
            cs = ax.pcolormesh(dlon,dlat,data,shading='nearest',
                        transform=ccrs.PlateCarree(),cmap=cmap,norm=norm)
        else:
            cs = ax.pcolormesh(dlon,dlat,data,shading='nearest',
                        transform=ccrs.PlateCarree(),cmap=cmap)
        cbar = fig.colorbar(cs,ax=ax,orientation='horizontal',pad=0,shrink=0.6)
        if bkcolor:
            cbar.set_label(label=clabel,size=16,color='white')
            ax.set_title(title,fontsize=18,color='white')
            cbar.ax.xaxis.set_tick_params(color='white', labelcolor='white')
        else:
            cbar.set_label(label=clabel,size=16)
            ax.set_title(title,fontsize=18)
        plt.tight_layout()

    if outpath:
        plt.savefig(outpath,dpi=400)
    plt.show()

#%% make map on a created axis
def mk_map(ax):
    """Function crates a map on cartopy axis.
    
    Parameters
    ----------
    ax: matplotlib figure axis
    
    Returns
    -------
    none.
    ------
    written by Katelyn O'Dell
    """

    ax.patch.set_visible(False)
    # plot shapfile with colors
    ax.add_feature(cfeature.LAND.with_scale('50m'),facecolor='gray',alpha=0.5)
    ax.add_feature(cfeature.OCEAN.with_scale('50m'))
    ax.add_feature(cfeature.STATES.with_scale('50m'),edgecolor='lightgray')
    ax.outline_patch.set_edgecolor('white')

#%% calculate wet pm2.5 at 35% RH and STP from wrf chem 
def calc_pmw(nc_fid,p25_str):
    """Function calculate wet pm2.5 at 35% RH and STP from wrf chem files from 
    Jian for GEO-XO project. Calculation based on email discussion between Jian He, 
    Brian McDonald,and Daven Henze and summarized in readme files 
    with the data from Jian.
    
    Parameters
    ----------
    nc_fid: loaded netCDF file ID
    
    Returns
    -------
    wet PM2.5 mass in ug/m3
    
    ------
    written by Katelyn O'Dell
    """
    so4a = nc_fid['so4a'][:]
    nh4a = nc_fid['nh4a'][:]
    no3a = nc_fid['no3a'][:]
    ec = nc_fid['ec'][:]
    orgpa = nc_fid['orgpa'][:]
    soa = nc_fid['soa'][:]
    p25 = nc_fid[p25_str][:]
    naa = nc_fid['naa'][:]
    cla = nc_fid['cla'][:]
    p = nc_fid['pres'][:]
    t = nc_fid['temp'][:]
    pm25w = (1.1*(so4a + nh4a + no3a) + ec + orgpa + 1.05*soa + p25 +1.86*(naa + cla))*(101325/p)*(t/298.0)
    return pm25w

#%% haversine function from Dr. Will Lassman
def haversine(lon0,lon1,lat0,lat1):
    """Function calculates distance between two lat/lon points assuming
    earth radius of 6,371 km
    
    Parameters
    ----------
    lon0: longitude of point 0 in degrees
    lon1: longitude of point 1 in degrees
    lat0: latitude of point 0 in degrees
    lat1: latitude of point 1 in degrees
    
    Returns
    -------
    Distance between point 0 and point 1 in meters.
    
    ------
    written by Will Lassman
    """

    r = 6371000. #m # mean earth radius                                                                                                                                                                                                                                              
    lon0 = lon0*np.pi/180

    lon1 = lon1*np.pi/180

    lat0 = lat0*np.pi/180

    lat1 = lat1*np.pi/180

    return 2*r*np.arcsin(np.sqrt(np.sin((lat1 - lat0)/2.)**2 +\
		 np.cos(lat0)*np.cos(lat1)*np.sin((lon1 - lon0)/2.)**2))

#%% acute HIA function
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
   
#%% chronic HIA function - generic
def chronic_HIA(conc, cf, pop, base_rate, betas, grid_area): 
    # beta is array of beta calculated from [rr,rr_lci,rr_uci]
    """Function description coming.
    
    Parameters
    ----------
    coming soon
    
    Returns
    -------
    coming soon
    
    ------
    written by Katelyn O'Dell
    """
    paf_avg_out = []
    events_tot_out = []
    events_tot_pp_out = []
    events_tot_pk_out = []
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

#%% gemm hia function       
def gemm_HIA(disease, theta, se_theta, alpha, mu, pi, base,
              bauPM, cvdPM,population):
    """Function description coming.
    
    Parameters
    ----------
    coming soon
    
    Returns
    -------
    coming soon
    
    ------
    written by Katelyn O'Dell, with significant help from Dr. Kelsey Bilsback
    """

    #print('Calc mortalities for ', disease)
    thetas = [theta - 2 * se_theta, theta, theta + 2 * se_theta]

    z_bau = np.where(bauPM > 2.4, bauPM - 2.4, 0)
    z_cvd = np.where(cvdPM > 2.4, cvdPM - 2.4, 0)
    
    Gamma_bau = np.log(1 + (z_bau / alpha)) / (1 + np.exp((mu - z_bau) / (pi)))
    Gamma_cvd = np.log(1 + (z_cvd / alpha)) / (1 + np.exp((mu - z_cvd) / (pi)))

    # Calculate hazard ratio
    HR_bau = np.exp(np.array(thetas)[:,None,None] * Gamma_bau)
    HR_cvd = np.exp(np.array(thetas)[:,None,None] * Gamma_cvd)

    #Mortalities
    M_bau = base * population * (1 - (1 / HR_bau))
    M_cvd = base * population * (1 - (1 / HR_cvd))

    dM_cvd = M_bau - M_cvd

    return(dM_cvd, M_bau, M_cvd)

#%% geo_idx from Dr. Gaige Kerr, https://github.com/gaigekerr/edf/blob/main/harmonize_afacs.py
def geo_idx(dd, dd_array):
    """Function searches for nearest decimal degree in an array of decimal 
    degrees and returns the index. np.argmin returns the indices of minimum 
    value along an axis. So subtract dd from all values in dd_array, take 
    absolute value and find index of minimum. 
    
    Parameters
    ----------
    dd : int
        Latitude or longitude whose index in dd_array is being sought
    dd_array : numpy.ndarray 
        1D array of latitude or longitude 
    
    Returns
    -------
    geo_idx : int
        Index of latitude or longitude in dd_array that is closest in value to 
        dd
    ------
    written by Dr. Gaige Kerr
    """
    geo_idx = (np.abs(dd_array - dd)).argmin()
    return geo_idx
#%%
def idwr(x, y, z, xi, yi):
    """Inverse distance weighting for interpolating gridded fields to census 
    tracts that are too small to intersect with the grid. 
    Parameters
    ----------
    x : list
        Longitude
    y : list
        Latitudes
    z : list
        NO2 column densities
    xi : list
        Unknown longitude
    yi : TYPE
        Unknown latitude
    Returns
    -------
    lstxyzi : list
        List comprised of unknown longitude, latitude, and interpolated 
        value for small census tract 
    ------
    written by Dr. Gaige Kerr, modified by Kate O'Dell to work with Will's haversine function
    """
    lstxyzi = []
    for p in range(len(xi)):
        lstdist = []
        for s in range(len(x)):
            # haversine inputs edited by Kate to match Will's haversine function
            d = (haversine(x[s], xi[p], y[s], yi[p]))
            lstdist.append(d)
        sumsup = list((1 / np.power(lstdist, 2)))
        suminf = np.nansum(sumsup)
        # The original configuration of this function had the following 
        # line of code
        # sumsup = np.nansum(np.array(sumsup) * np.array(z))
        # However, there were issues with arrays with missing data and NaNs, 
        # so it was changed to the following: 
        sumsup = np.nansum(np.array(sumsup) * np.array(z))
        u = sumsup / suminf
        xyzi = [xi[p], yi[p], u]
        lstxyzi.append(xyzi)
    return lstxyzi
    