#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 20 16:22:58 2021

Global intrinsic unc time series 

@author: ccamargo
"""

import numpy as np
import xarray as xr
import pandas as pd
import sys
sys.path.append("/Users/ccamargo/Documents/py_scripts/")
import utils_SL as sl 
import utils_hec as hec
import os
# import cmocean as cmo
# from numba import jit
import datetime as dt

import utils_SLE_v2 as sle
import matplotlib.pyplot as plt
#%% GRACE
path= '/Volumes/LaCie_NIOZ/data/barystatic/intrinsic_unc/use/'
grace = 'GRA-Mascon-JPL_nobuf_sel_180x360.nc'

ds=xr.open_dataset(path+grace)
ds=ds.drop_sel(name=['TCWS_JPL'])
#% % namelist
print(ds.name)
names= [name for name in np.array(ds.name) ]
time1=np.array(ds.time)
tdec,_=sl.get_dec_time(time1)
names2=[]
for name in names:
    if name.startswith('GLWS'):
        name = name.replace('GLWS','GLA')
    if len(name.split('_'))==3:
        if name.endswith('gl'):
            name=name.split('_')[0]+'_'+name.split('_')[1]
        else:
            name = name.split('_')[0]+'_'+name.split('_')[2]
    print(name)
    names2.append(name)
        
names
len_reg=len(names2)
# #%% regional to global
# glb_mm = np.zeros((len(ds.name),len(ds.time)))
# # glb_EWH = np.full_like(glb_mm,0)

# # glb_mm_y = np.zeros((len(ds.name),len(ds.year)))
# # # glb_EWH_y = np.full_like(glb_mm_y,0)
# lat=np.array(ds.lat)
# lon=np.array(ds.lon)
# for iname in range(0, len(ds.name)):
#     for itime in range(0,len(ds.time)):
#         glb_mm[iname,itime]=sle.reg_to_glb(ds.SL_mm[iname,itime,:,:], lat, lon)[0]
#%% namelist
print(ds.name)
names= [name for name in np.array(ds.name) ]

names2=[]
for name in names:
    if name.startswith('GLWS'):
        name = name.replace('GLWS','GLA')
    if len(name.split('_'))==3:
        if name.endswith('gl'):
            name=name.split('_')[0]+'_'+name.split('_')[1]
        else:
            name = name.split('_')[0]+'_'+name.split('_')[2]
    print(name)
    names2.append(name)
        
names
len_reg=len(names2)

ds['name']=names2

#%% compute ocean mean

df_glb=pd.DataFrame(np.hstack(np.array(tdec)),columns=['time'])
for iname, name in enumerate(names2):
    ds3 = ds.sel(name=name)
    
    value=np.array(ds3.SL_mm[:,:,:])
    mu_ts = - np.nansum(value.reshape(value.shape[0],value.shape[1]*value.shape[2]),axis=1)
    # we put the negative value because mass loss = SL gain
    
    # replace 0 (missing values) for NaN:
    mu_ts[np.where(mu_ts==0)]=np.nan
    # compute the mean over time:
    mu=np.nanmean(mu_ts)
    df_glb[name]=mu_ts-mu
    df_glb[name].plot();
    plt.title(name);
    # plt.show()
    plt.close()
#%% IMBIE
path1 ='/Volumes/LaCie_NIOZ/data/barystatic/source_var/'
file2 = 'means_v2.nc'
ds=xr.open_dataset(path1+file2)
# data=np.array(ds.GIS_UCI_slc)
# plt.plot(data[:,-2])

ds = ds.sel(time=slice(1993,ds.time[-1]))

time2=np.array(ds.time)

#% % aadd to namelist
names = ['GIS_IMB','AIS_IMB']
glb=np.zeros((len(ds.time),2))
j=0
jj=0
for iname,name in enumerate(names):
    da=ds[name+'_unc']
    # print(name)
    if name.split('_')[0] in da['reg_'+name]:
        if name.split('_')[-1]=='IMB':
            glb[:,j]=np.array(da.data[0,:])
            j=j+1
    

#%% combine with regional:
names_tot=[name for name in names2]
names_tot.extend(names)
glb_mm_month=np.zeros((len(time1),len(names_tot)))
# glb_mm_month[0:len_reg,:]=glb_ts_mm
glb_mm_month[:,0:len_reg]=np.array(df_glb[names2])
glb_mm_month[0:len(time2),-2]=glb[:,0]
glb_mm_month[0:len(time2),-1]=glb[:,1]

#%%
glb_mm_month[np.where(glb_mm_month==0)]=np.nan

#%%
for i,n in enumerate(names_tot):
    plt.plot(glb_mm_month[:,i],label=n)
    plt.legend()
    
    plt.close()



#%% 
da = xr.Dataset(data_vars={'OM_monthly_intrinsic_unc':(('months','name_monthly'),glb_mm_month),
                           
                           },
                         coords={'name_monthly':names_tot,
                                 'months':time1,
                                 })

da.attrs['units']='mm'
da.attrs['Metadata']='Ocean mean intrinsic unc SLC variations from different sources (AIS,GIS,GLA,LWS) and different datasets'
da.attrs['script']='7d-OM_global_intrinsic.py'

path='/Volumes/LaCie_NIOZ/data/barystatic/global/'
da.to_netcdf(path+'global_intrinsic-timeseries_v2.nc')


