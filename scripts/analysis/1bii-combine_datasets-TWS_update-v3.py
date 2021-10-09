#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 16:19:01 2021

@author: ccamargo
"""

import numpy as np
# import scipy.optimize as opti
import xarray as xr
import matplotlib.pyplot as plt
import sys
sys.path.append("/Users/ccamargo/Documents/py_scripts/")
import utils_SL as sl 
import utils_SLE_v2 as sle

# from netCDF4 import Dataset
import pandas as pd
import os

import datetime as dt

import cmocean as cm
# from mpl_toolkits.basemap import Basemap
# from matplotlib.gridspec import GridSpec
from cartopy import crs as ccrs#, feature as cfeature

#%%  Open dataset made in barystatic_standardize2.py
path='/Volumes/LaCie_NIOZ/data/barystatic/use/'
ds=xr.open_dataset(path+'comb/ALL_datasets_1993-2020_180x360_v2.nc')

# We need to correct PCR-GLOBWB
# the unit is METERS of Equivalent Water Thickness [x1000]
# and not dcm of Equivalent Water Thickness [x100]
file='/Volumes/LaCie_NIOZ/data/barystatic/use/TWS_PCR-GLOBWB_vW5E5_noGLA_180x360.nc'
name='LWS_GWB'
SL_mm=np.array(ds.SL_mm)
SL_mm_y=np.array(ds.SL_mm_y)
SL_EWH=np.array(ds.SL_EWH)
SL_EWH_y=np.array(ds.SL_EWH_y)

    # SL_mm     (name, time, lat, lon) float64 ...
    # SL_EWH    (name, time, lat, lon) float64 ...
    # SL_mm_y   (name, year, lat, lon) float64 ...
    # SL_EWH_y  (name, year, lat, lon) float64 ...
    

#%% time vectors
# time range: months 1993-2020.08
time,time2=sl.make_date('01-01-1993',(27*12)+8)
tdec,tdec0=sl.get_dec_time(time2,ns=1e-6)
t=[dt.datetime.utcfromtimestamp(t.astype(int) * 1e-6) for t in time2]
ty=np.array([dt.datetime.utcfromtimestamp(t.astype(int) * 1e-6).timetuple().tm_year
                                for t in time2])
tm=np.array([dt.datetime.utcfromtimestamp(t.astype(int) * 1e-6).timetuple().tm_mon
                                for t in time2])

# yearly
timey=np.arange(1993,2021)
idx=np.zeros((len(timey)))
for i,year in enumerate(timey):
    idx[i]=np.where(ty==year)[0][0]
    
#%% TWS GWB
ifile = 24
iname=ifile
print(ds.name[ifile])
if name==ds.name[ifile]:
    da=xr.open_dataset(file)
    da_y=da.groupby('time.year').mean(dim='time')
    name='tws_no_gla'

    #% % find time index
    time_local=np.array(da.time)
    tdec_local,tdec0=sl.get_dec_time(time_local,ns=1e-9)
    lon=np.array(da.lon)
    lat=np.array(da.lat)
    timey_local=np.array(da_y.year)
    
    # find index for yearly:
    idx_local=np.zeros((len(timey)))
    idx_local.fill('nan')
    # idx_local=np.full_like(timey,np.nan)
    for i,year in enumerate(timey):
        if year<=np.max(timey_local) and \
            year>=np.min(timey_local):
            idx_local[i]=np.where(timey_local==year)[0][0]
    
    idx_0y=np.array(np.where(np.isfinite(idx_local)))[0,0]
    idx_00y=int(idx_local[np.isfinite(idx_local)][0])
    idx_1y=np.array(np.where(np.isfinite(idx_local)))[0,-1]+1
    
    
    # find index for monthly:
    ty_local=np.array([dt.datetime.utcfromtimestamp(t.astype(int) * 1e-9).timetuple().tm_year
                                    for t in time_local])
    tm_local=np.array([dt.datetime.utcfromtimestamp(t.astype(int) * 1e-9).timetuple().tm_mon
                                    for t in time_local])
    
    idx_local_m=np.zeros((len(ty)))
    idx_local_m.fill('nan')
    t2_local=[dt.datetime.utcfromtimestamp(t.astype(int) * 1e-9 ) for t in time_local]
    for i,iy in enumerate(t):
        year=iy.year
        if np.any(ty_local==year):
            month=iy.month
            if np.any(tm_local==month):
                # print(year)
                # print(month)
                jdx=np.array(np.where(year==ty_local))[0]
                jdx2=np.array(np.where(month==tm_local))[0]
                for jd in jdx2:
                    if np.any(jdx==jd):
                        j=jdx[jdx==jd]
                        idx_local_m[i]=j
    
    #% % Water thickness original (m of water thickness)
    data=np.array(da[name])*1000 # m->mm of water thickness
    datay=np.array(da_y[name])*1000 
    # Put data in main array
    SL_EWH_y[iname,idx_0y:idx_1y,:,:]=datay[idx_00y:len(timey_local),:,:]
    for i,j in enumerate(idx_local_m):
        if np.isfinite(j):
            j=int(j)
            SL_EWH[iname,i,:,:]=data[j,:,:]
    # Transform data: 
    for i in range(len(data)):
        data[i,:,:]=sle.EWH_to_height(data[i,:,:]).reshape(180,360)
    for i in range(len(datay)):
        datay[i,:,:]=sle.EWH_to_height(datay[i,:,:]).reshape(180,360)
    # Put transformed data in main array:
    SL_mm_y[iname,idx_0y:idx_1y,:,:]=datay[idx_00y:len(timey_local),:,:]
    for i,j in enumerate(idx_local_m):
        if np.isfinite(j):
            j=int(j)
            SL_mm[iname,i,:,:]=data[j,:,:]

#%% make and save data array
ds_mask=xr.open_dataset('/Volumes/LaCie_NIOZ/data/barystatic/masks/barystatic_mask2.nc')
ds_mask

#%%
da=xr.Dataset(data_vars={'SL_mm':(('name','time','lat','lon'),SL_mm),
                         'SL_EWH':(('name','time','lat','lon'),SL_EWH),
                         'SL_mm_y':(('name','year','lat','lon'),SL_mm_y),
                         'SL_EWH_y':(('name','year','lat','lon'),SL_EWH_y),
                         'mask':(('lat','lon'),ds_mask['mask6']),
                         'mask2':(('lat','lon'),ds_mask['mask12']),
                         
                        },
                         coords={'lat':ds.lat,
                                 'lon':ds.lon,
                                 'time':time2,
                                 'tdec':tdec,
                                 'year':timey,
                                 'name':ds.name})
da['SL_mm'].attrs['units']='mm of sea level'
da['SL_mm_y'].attrs['units']='mm of sea level'
da['SL_EWH'].attrs['units']='mm of Equivalent Water Thickness'
da['SL_EWH_y'].attrs['units']='mm of Equivalent Water Thickness'
da['SL_mm'].attrs['long_name']='Monthly ocean mass change in mm of sea level height'
da['SL_EWH'].attrs['long_name']='Monthly ocean mass change in mm of equivalent water thickness'
da['SL_mm_y'].attrs['long_name']='Yearly averages of ocean mass change in mm of sea level height'
da['SL_EWH_y'].attrs['long_name']='Yearly averages of ocean mass change in mm of equivalent water thickness'

da.attrs=ds.attrs
da.attrs['script']='barystatic_standardize_update.py'
da.attrs['date_created']=str(dt.datetime.now())
da.to_netcdf(path+'comb/ALL_datasets_1993-2020_180x360_v3_update.nc')

#%% Selection updated
path='/Volumes/LaCie_NIOZ/data/barystatic/use/comb/'
da=xr.open_dataset(path+'ALL_datasets_1993-2020_180x360_v3_update.nc')
print(da.name)

da=da.sel(name=[#'AIS_IMB', 'AIS_R19_basins',
                'GLWS_ZMP', 'GLWS_WGP_gl', 
                'AIS_300_CSR', 'GIS_300_CSR', 'LWS_CSR', 'GLWS_CSR', 
                # 'TCWS_CSR', 
                # 'AIS_CSR', 'GIS_CSR',
                'AIS_300_JPL', 'GIS_300_JPL', 'LWS_JPL', 'GLWS_JPL', 
                # 'TCWS_JPL',
                # 'AIS_JPL', 'GIS_JPL', 
                # 'AIS_proj_CSR', 'GIS_proj_CSR', 'AIS_proj_JPL','GIS_proj_JPL', 
                # 'GIS_IMB', 'GIS_M19', 
                'LWS_GWB', 'LWS_WGP_gl',
                
               # 'TCWS_WaterGAP'
               ])
print(da)
da.to_netcdf(path+'12_input_datasets_buf_1993-2020_180x360_update_v3.nc')

#%%
#%% Selection updated
path='/Volumes/LaCie_NIOZ/data/barystatic/use/comb/'
da=xr.open_dataset(path+'ALL_datasets_1993-2020_180x360_v3_update.nc')
print(da.name)

da=da.sel(name=['AIS_IMB', 'AIS_R19_basins',
                'GLWS_ZMP', 'GLWS_WGP_gl', 
                'AIS_300_CSR', 'GIS_300_CSR', 'LWS_CSR', 'GLWS_CSR', 
                # 'TCWS_CSR', 
                # 'AIS_CSR', 'GIS_CSR',
                'AIS_300_JPL', 'GIS_300_JPL', 'LWS_JPL', 'GLWS_JPL', 
                # 'TCWS_JPL',
                # 'AIS_JPL', 'GIS_JPL', 
                # 'AIS_proj_CSR', 'GIS_proj_CSR', 'AIS_proj_JPL','GIS_proj_JPL', 
                 'GIS_IMB', 'GIS_M19', 
                'LWS_GWB', 'LWS_WGP_gl',
                
               # 'TCWS_WaterGAP'
               ])
print(da)
da.to_netcdf(path+'16_input_datasets_buf_1993-2020_180x360_update_v3.nc')