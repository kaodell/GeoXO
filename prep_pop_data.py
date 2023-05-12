# -*- coding: utf-8 -*-
"""
prep_pop_data.py
    python script to regrid population data to the WRF Chem grid
    adapted from https://github.com/kaodell/smoke_HIA/blob/main/regrid_population.py
    02.03.22 - modified to regrid to Shobha's 4 km grid
    03.03.22 - modfied back to the 12 km grid
    03.08.22 - modified to make 2020 pop data for Molly Robertson on the kriging grid
    04.14.23 - final version for GeoXO project
"""

#%% user inputs
# indicate file locations
# project folder
prj_folder = '/Users/kodell/Library/CloudStorage/GoogleDrive-kodell@email.gwu.edu/My Drive/Ongoing Projects/GeoXO/'
# grid to regrid population to
to_grid_fn = '/Users/kodell/Library/CloudStorage/Box-Box/Shobha_data/new_ABI/ABI_daily/2020/ABI_CONUS_PM2.5_daily_20200101_0.01x0.01.nc'
# area of the to-grid
to_grid_area_fn = 'concentration_data/ACX_analysis/abi_corner_area.nc'
# population density file to regrid
pop_fn = 'population_data/NASA_SEDAC_pop/gpw-v4-population-density-rev11_totpop_2pt5_min_nc/gpw_v4_population_density_rev11_2pt5_min.nc'

# US national shapefile to sum US pop to check we are close and to make masks for plotting
shp_fn = 'population_data/cb_2018_us_nation_5m/cb_2018_us_nation_5m'

# US states shapefile to make a california mask for calculating california totals
cal_shp_fn = prj_folder + 'population_data/tl_2020_us_state/tl_2020_us_state'

# name of output file
out_fn = 'population_data/NASA_SEDAC_pop/regridded/regrided_2020pop_0.01_final.nc'

#%% load modules
import numpy as np
import netCDF4
import datetime
from mpl_toolkits.basemap import interp
import shapefile
from ODell_udf import plt_map
import matplotlib as mplt

#%% load 'to' grid
nc_fid = netCDF4.Dataset(to_grid_fn)
to_lat = nc_fid.variables['lat'][:]
to_lon = nc_fid.variables['lon'][:]
nc_fid.close()

# replace fill values with nans
to_lat = np.where(to_lat.data==to_lat.fill_value,np.nan,to_lat.data)
to_lon = np.where(to_lon.data==to_lon.fill_value,np.nan,to_lon.data)

# create meshgrid
to_glon, to_glat = np.meshgrid(to_lon,to_lat)

# to grid area
# calc area
height_m = 1111.1 #m height is the same everywhere
width_deg = 0.01
widths_m = np.empty(to_glon.shape)
for j in range(len(to_lat)):
    lat = to_lat[j]*np.pi/180.0 # convert to radians
    widths_m[j,:] = np.cos(lat) * 1111.1 # distance of 0.01 degree in meters at equator
grid_area = (height_m/1000.0)*(widths_m/1000.0) # convert to km

#%% load population data and 'from' grid
# population density from SEDAC
nc_fid = netCDF4.Dataset(prj_folder+pop_fn) # open file to read

# data to regrid
pop_density_data = nc_fid.variables['Population Density, v4.11 (2000, 2005, 2010, 2015, 2020): 2.5 arc-minutes'][:] # people per km2 
og_lon = nc_fid.variables['longitude'][:].data
og_lat = nc_fid.variables['latitude'][:].data
nc_fid.close()

pop_density_wm = pop_density_data[4] # 4 = 2020
land_area_wm = pop_density_data[8] + pop_density_data[9] # add water for pop counts to match the regridded calcs
national_ids = pop_density_data[10]
# replace fill value with zero
pop_density = np.ma.getdata(pop_density_wm)
land_area = np.ma.getdata(land_area_wm)
mask = np.ma.getmask(pop_density_wm)
la_mask = np.ma.getmask(land_area_wm)
pop_density[mask] = 0.0
land_area[mask] = np.nan

# cut to the US
US_ind_latmax = np.max(np.where(og_lat > 20))
US_ind_latmin = np.min(np.where(og_lat < 60))
US_ind_lonmax = np.max(np.where(og_lon < -60))
US_ind_lonmin = np.min(np.where(og_lon > -130))

US_og_lon = np.copy(og_lon[US_ind_lonmin:US_ind_lonmax])
US_og_lat = np.copy(og_lat[US_ind_latmin:US_ind_latmax])
US_pop_density = np.copy(pop_density[ US_ind_latmin:US_ind_latmax,US_ind_lonmin:US_ind_lonmax])
US_land_area = np.copy(land_area[ US_ind_latmin:US_ind_latmax,US_ind_lonmin:US_ind_lonmax])
US_ids = np.copy(national_ids[ US_ind_latmin:US_ind_latmax,US_ind_lonmin:US_ind_lonmax])

# plot grid areas to check
US_og_lon_m, US_og_lat_m = np.meshgrid(US_og_lon,US_og_lat)
plt_map(US_og_lon_m,US_og_lat_m,US_land_area,1,'jet','grid area [km]','from grid area',clim=[10,20])
plt_map(to_glon,to_glat,grid_area,1,'jet','grid area [km]','to grid area check')

#%% run grid interpolation
# lat has to be increasing, so flip
US_og_lat_flip = np.flipud(US_og_lat)
US_pop_density_flip = np.flipud(US_pop_density)
pop_density_regrid = interp(US_pop_density_flip, US_og_lon, US_og_lat_flip, to_glon, to_glat, order=0) 
# order=0 is a nearest neighbor interpolation

# use grid area above to get counts
pop_regrid = pop_density_regrid * grid_area

#%% make figures to compare
# population density
plt_map(US_og_lon_m, US_og_lat_m, US_pop_density,0.1, 'magma', 'population density', 
        'NASA SEDAC pop density 2020',clim=[0,100])
plt_map(to_glon, to_glat, pop_density_regrid,0.1, 'magma', 'population density', 
        'regridded pop density 2020',clim=[0,100])

#%% print US pop totals with both grids to check
# first we need a US mask for the 'to' grid
# load shapefile
shps_file = shapefile.Reader(prj_folder+shp_fn)
shps_shp = shps_file.shape(0)
shps_records = shps_file.records()
shps_shapes = shps_file.shapes()
# loop through shapes and create mask
si = 0
area_mask = np.zeros(to_glon.shape)
for j in range(len(shps_records)):
        area_shp = shps_shapes[j]
        for i in range(len(area_shp.parts)):
            i0 = area_shp.parts[i]
            if i < len(area_shp.parts)-1:
            		i1 = area_shp.parts[i+1] - 1
            else:
            		i1 = len(area_shp.points)
            seg = area_shp.points[i0:i1+1]
            mpath = mplt.path.Path(seg)
            points = np.array((to_glon.flatten(), to_glat.flatten())).T
            mask = mpath.contains_points(points).reshape(to_glon.shape)
            area_inds = np.where(mask==True)
            area_mask[area_inds] = 1
# check mask
plt_map(to_glon, to_glat, area_mask,1, 'jet', 'mask', 'US border mask',
        clim=[0,1])
# pull US inds from the 'from' grid with the country IDs
us_inds = np.where(US_ids == 840)
us_pop_og = US_land_area[us_inds[0],us_inds[1]]*US_pop_density[us_inds[0],us_inds[1]]

# total pop check
print('estimated contig. us pop regrid',np.sum(area_mask*pop_regrid))           
print('estimated contig. us pop SEDAC',np.sum(us_pop_og))           
print(100.0*(np.sum(area_mask*pop_regrid)-np.sum(us_pop_og))/np.sum(us_pop_og),'% difference')

#%% also make california mask for calculating california totals
# load California shape file
shps_file = shapefile.Reader(cal_shp_fn)
shps_shp = shps_file.shape(0)
shps_records = shps_file.records()
shps_shapes = shps_file.shapes()
# find california shape
si = 0
cal_area_mask = np.zeros(to_glon.shape)
for j in range(len(shps_records)):
    if shps_records[j][2] == '06':
        print(shps_records[j][2])
        area_shp = shps_shapes[j]
for i in range(len(area_shp.parts)):
    i0 = area_shp.parts[i]
    if i < len(area_shp.parts)-1:
    		i1 = area_shp.parts[i+1] - 1
    else:
    		i1 = len(area_shp.points)
    seg = area_shp.points[i0:i1+1]
    mpath = mplt.path.Path(seg)
    points = np.array((to_glon.flatten(), to_glat.flatten())).T
    mask = mpath.contains_points(points).reshape(to_glon.shape)
    area_inds = np.where(mask==True)
    cal_area_mask[area_inds] = 1
        
# check mask
plt_map(to_glon, to_glat, cal_area_mask,1, 'jet', 'mask', 'cal border mask',
        clim=[0,1])

#%% write to file
# create netCDF file
nc_w_fid = netCDF4.Dataset(prj_folder+out_fn, 'w', format='NETCDF4')
nc_w_fid.description = 'Population from SEDAC in 2020'
nc_w_fid.history = 'Created' + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

# define file dimensions
nc_w_fid.createDimension('time', None) # unlimited dimension
nc_w_fid.createDimension('gridx', to_glon.shape[0])
nc_w_fid.createDimension('gridy', to_glon.shape[1])

lat_w = nc_w_fid.createVariable('lat', np.float32, ('gridx',))
lon_w = nc_w_fid.createVariable('lon', np.float32, ('gridy',))
pop_w = nc_w_fid.createVariable('population', np.float32, ('gridx','gridy',))
pop_density_w = nc_w_fid.createVariable('population_density', np.float32, ('gridx','gridy',))
grid_area_w = nc_w_fid.createVariable('grid_area', np.float32, ('gridx','gridy',))
us_mask_w = nc_w_fid.createVariable('us_mask', np.float32, ('gridx','gridy',))
cal_mask_w = nc_w_fid.createVariable('cal_mask', np.float32, ('gridx','gridy',))

lon_w[:] = to_lon
lat_w[:] = to_lat
pop_w[:] = pop_regrid
pop_density_w[:] = pop_density_regrid
grid_area_w[:] = grid_area
us_mask_w[:] = area_mask
cal_mask_w[:] = cal_area_mask

nc_w_fid.close()

