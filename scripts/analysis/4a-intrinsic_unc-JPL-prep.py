#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  7 11:07:47 2021

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


#%% open mask
ds=xr.open_dataset('/Volumes/LaCie_NIOZ/data/barystatic/masks/barystatic_mask2.nc')
ds
ds.gre.plot()
ds.mask3.plot()

mask=np.array(ds.mask3)
# '0=ocean,1=land,2=gre 300km buffer ,3=ant 300km buf,4=glaciers'
regions=['ocean','tws','gre','ant','gla','bary_total']
# maskant=np.copy(mask)
# maskant[np.where(mask==3)]=100
# maskant[np.where(maskant<100)]='nan'

#%% mascons
# dataset=['JPL']
# # path='/Volumes/LaCie_NIOZ/data/barystatic/original/GRA-Mascon-'
# path_to_save='/Volumes/LaCie_NIOZ/data/barystatic/intrinsic_unc/GRA-Mascon-'
# path='/Volumes/LaCie_NIOZ/data/barystatic/regrid/GRA-Mascon-'
# for d in dataset:
#     fin=path+d+'_180x360.nc'
#     ds=xr.open_dataset(fin)
#     print(ds)
    
#     time=np.arange(0,len(ds.time))

    
#     for i,reg in  enumerate(regions):
#         data=np.zeros((len(time),len(ds.lat),len(ds.lon)))
#         # if d=='JPL':
#         #     data_unc=np.zeros((len(time),len(ds.lat),len(ds.lon)))
            
#         for j in time:
#             tmp=np.array(ds.uncertainty[j,:,:])
#             if i<5:
#                 tmp[np.where(mask!=i)]=np.nan
#             else:
#                 tmp[np.where(mask==0)]=np.nan
#             data[j,:,:]=tmp
            
#             # if d=='JPL':
#             #     tmp_unc=np.array(ds.uncertainty[j,:,:])
#             #     tmp_unc[np.where(mask!=i)]='nan'
#             #     data_unc[j,:,:]=tmp_unc
        
#         # out 
#         ds[reg]=(('time','lat','lon'),data)
#         # if d=='JPL':
#         #     ds[reg+'_unc']=(('time','lat','lon'),data_unc)    
    
#     ds['mask']=(('lat','lon'),mask)
#     ds.attrs['regions']='Regions selected from the mascons using the barystatic component mask'
#     ds.attrs['bary_mask']='0=ocean,1=land,2=gre 300km buffer ,3=ant 300km buf,4=glaciers'
#     ds.attrs['script']='isolate_components_mascon.py'
#     ds.attrs['missing_value']='masked areas were replaced by 9999'
#     ds.to_netcdf(path_to_save+d+'_300kmbuf_sel_180x360.nc')
    

#%% without buffer
ds=xr.open_dataset('/Volumes/LaCie_NIOZ/data/barystatic/masks/barystatic_mask2.nc')
ds
# ds.gre.plot()
# ds.mask3.plot()

mask=np.array(ds.mask2)
# '0=ocean,1=land,2=gre 300km buffer ,3=ant 300km buf,4=glaciers'
regions=['ocean','tws','gre','ant','gla','bary_total']
# maskant=np.copy(mask)
# maskant[np.where(mask==3)]=100
# maskant[np.where(maskant<100)]='nan'

#%% mascons
dataset=['JPL']
# path='/Volumes/LaCie_NIOZ/data/barystatic/original/GRA-Mascon-'
path_to_save='/Volumes/LaCie_NIOZ/data/barystatic/intrinsic_unc/GRA-Mascon-'
path='/Volumes/LaCie_NIOZ/data/barystatic/regrid/GRA-Mascon-'
for d in dataset:
    fin=path+d+'_180x360.nc'
    ds=xr.open_dataset(fin)
    print(ds)
    
    time=np.arange(0,len(ds.time))

    
    for i,reg in  enumerate(regions):
        data=np.zeros((len(time),len(ds.lat),len(ds.lon)))
        # if d=='JPL':
        #     data_unc=np.zeros((len(time),len(ds.lat),len(ds.lon)))
            
        for j in time:
            tmp=np.array(ds.uncertainty[j,:,:])
            if i<5:
                tmp[np.where(mask!=i)]=np.nan
            else:
                tmp[np.where(mask==0)]=np.nan
            data[j,:,:]=tmp
            
            # if d=='JPL':
            #     tmp_unc=np.array(ds.uncertainty[j,:,:])
            #     tmp_unc[np.where(mask!=i)]='nan'
            #     data_unc[j,:,:]=tmp_unc
        
        # out 
        ds[reg]=(('time','lat','lon'),data)
        # if d=='JPL':
        #     ds[reg+'_unc']=(('time','lat','lon'),data_unc)    
    
    ds['mask']=(('lat','lon'),mask)
    ds.attrs['regions']='Regions selected from the mascons using the barystatic component mask'
    ds.attrs['bary_mask']='0=ocean,1=land,2=gre 300km buffer ,3=ant 300km buf,4=glaciers'
    ds.attrs['script']='isolate_components_mascon.py'
    ds.attrs['missing_value']='masked areas were replaced by nans'
    ds.to_netcdf(path_to_save+d+'_nobuf_sel_180x360.nc')
    


# %% standardize time

flist=sl.get_filelist(path_to_save,ext='*.nc')
# print(flist)
# flist=['/Volumes/LaCie_NIOZ/data/barystatic/use/GRA-Mascon-CSR_noBuf_sel_180x360.nc',
#       '/Volumes/LaCie_NIOZ/data/barystatic/use/GRA-Mascon-JPL_noBuf_sel_180x360.nc']
names=[
        'AIS_JPL','GIS_JPL',
        'LWS_JPL','GLWS_JPL',
        'TCWS_JPL',  
      ] 

variables = ['ant','gre','tws','gla','bary_total']

ds_mask=xr.open_dataset('/Volumes/LaCie_NIOZ/data/barystatic/masks/barystatic_mask2.nc')

#% % make dates
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
    


#% %
for f in flist:
    ds=xr.open_dataset(f)
    ds_y=ds.groupby('time.year').mean(dim='time')
    
    #% % make empty arrays:
    SL_mm=np.zeros((len(names),len(time),180,360))
    SL_mm.fill('nan')
    SL_EWH=np.full_like(SL_mm,np.nan)
    
    SL_mm_y=np.zeros((len(names),len(timey),180,360))
    SL_mm_y.fill('nan')
    SL_EWH_y=np.full_like(SL_mm_y,np.nan)

    # % % find time index
    time_local=np.array(ds.time)
    tdec_local,tdec0=sl.get_dec_time(time_local,ns=1e-9)
    lon=np.array(ds.lon)
    lat=np.array(ds.lat)
    timey_local=np.array(ds_y.year)
    
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
    i=225
    idx_local_m[i+1]=idx_local_m[i]
    idx_local_m[i]=idx_local_m[i]-1
    i=267
    idx_local_m[i+1]=idx_local_m[i]
    idx_local_m[i]=idx_local_m[i]-1

    # #############
    # loop over variables 
    for ivar, var in enumerate(variables):
        # # Water thickness
        data=np.array(ds[var])*10 # mm of water thickness
        datay=np.array(ds_y[var])*10 
        # Put data in main array
        SL_EWH_y[ivar,idx_0y:idx_1y,:,:]=datay[idx_00y:len(timey_local),:,:]
        for i,j in enumerate(idx_local_m):
            if np.isfinite(j):
                j=int(j)
                SL_EWH[ivar,i,:,:]=data[j,:,:]
        # Transform data: 
        for i in range(len(data)):
            data[i,:,:]=sle.EWH_to_height(data[i,:,:]).reshape(180,360)
        for i in range(len(datay)):
            datay[i,:,:]=sle.EWH_to_height(datay[i,:,:]).reshape(180,360)
        # Put transformed data in main array:
        SL_mm_y[ivar,idx_0y:idx_1y,:,:]=datay[idx_00y:len(timey_local),:,:]
        for i,j in enumerate(idx_local_m):
            if np.isfinite(j):
                j=int(j)
                SL_mm[ivar,i,:,:]=data[j,:,:]
   
    #% % make and save data array
    
    #% %
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
                                      'name':names})
    da['SL_mm'].attrs['units']='mm of sea level'
    da['SL_mm_y'].attrs['units']='mm of sea level'
    da['SL_EWH'].attrs['units']='mm of Equivalent Water Thickness'
    da['SL_EWH_y'].attrs['units']='mm of Equivalent Water Thickness'
    da['SL_mm'].attrs['long_name']='Monthly ocean mass change in mm of sea level height'
    da['SL_EWH'].attrs['long_name']='Monthly ocean mass change in mm of equivalent water thickness'
    da['SL_mm_y'].attrs['long_name']='Yearly averages of ocean mass change in mm of sea level height'
    da['SL_EWH_y'].attrs['long_name']='Yearly averages of ocean mass change in mm of equivalent water thickness'
    # da['time'].attrs['long_name']='Time'
    # da['time'].attrs['standard_time']='time'
    # da['time'].attrs['units']='months since january-1993'
    # da['time'].attrs['calendar']='gregorian'
    # da['tdec'].attrs['long_name']='Decimal time'
    # da['tdec'].attrs['units']='months since january-1993'
    # da['tdec'].attrs['calendar']='gregorian'
    # da['year'].attrs['long_name']='Years'
    # da['year'].attrs['datasets']='R19,M19 and ZMP were originally yearly data. Other datasets have been averaged to obtain an year mean.'
    # da['tdec'].attrs['units']='years 1993'
    da['name'].attrs['long_name']='Name of the source region/contribution of ocean mass change'
    #% %
    da.attrs['metadata']='Ocean mass sea level changes (in mm and EWH), from differente contributions (ANT, GRE, TWS, GLA) \
          from JPL Mascons)'
    
    da.attrs['Author']='Carolina M.L. Camargo'
    da.attrs['reference']='Camargo et al. 2021'
    da.attrs['date_created']=str(dt.datetime.now())
    
    path='/Volumes/LaCie_NIOZ/data/barystatic/intrinsic_unc/use/'
    da.to_netcdf(path+f.split('/')[-1])
#%%

