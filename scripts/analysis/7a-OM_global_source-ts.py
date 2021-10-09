#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 23 15:43:35 2021

 Global mean time series -
 Computed by just summing up the land mass variations
 This was confirmed with the mean value printed out from the SLE
 
@author: ccamargo
"""



import numpy as np
import xarray as xr
import sys
sys.path.append("/Users/ccamargo/Documents/py_scripts/")
import utils_SL as sl 
import utils_hec as hec
import os
# import cmocean as cmo
# from numba import jit
import datetime as dt
import pandas as pd

import utils_SLE_v2 as sle
import matplotlib.pyplot as plt

#%% regional

path ='/Volumes/LaCie_NIOZ/data/barystatic/source_var/'
file1='regional_v3.nc'

ds=xr.open_dataset(path+file1)
ds=ds.drop_sel(name=['TCWS_CSR','TCWS_JPL'])
time1=np.array(ds.time)
years1=np.array(ds.year)
tdec,tdec0=sl.get_dec_time(time1)
# ds['time']=tdec
# t0=2003
# t1=2017
# ds= ds.sel(time=slice(t0,t1))
# tdec=np.array(ds.time)
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

lat=np.array(ds.lat)
lon=np.array(ds.lon)
dimlon=len(lon)
dimlat=len(lat)
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
    # out= sl.get_ts_trend(np.array(df_glb['time']),np.array(df_glb[name]),offset=1,plot=False)
    
    # print(name)
    # print(out[0])
    
    # out=sl.get_reg_trend(np.array(df_glb['time']),value, lat,lon)
    # plt.pcolor(out[0]);plt.title(name)
    # print(np.nansum(out[0]))
    
    # out=sle.run_SLE(sle.height_to_vol(out[0]).reshape(dimlat,dimlon),'_test')

#%% basin means
file2 = 'means_v2.nc'
ds=xr.open_dataset(path+file2)
# data=np.array(ds.GIS_UCI_slc)
# plt.plot(data[:,-2])

ds = ds.sel(time=slice(1993,ds.time[-1]))
ds = ds.sel(year=slice(1993,ds.year[-1]))

time2=np.array(ds.time)
years2=np.array(ds.year)


 #%% aadd to namelist
names = ['GIS_IMB','AIS_IMB','AIS_UCI','GIS_UCI','GLA_ZMP']
glb=np.zeros((len(ds.time),2))
glb_y=np.zeros((len(ds.year),3))
j=0
jj=0
for iname,name in enumerate(names):
    da=ds[name+'_slc']
    # print(name)
    if name.split('_')[0] in da['reg_'+name]:
        if name.split('_')[-1]=='IMB':
            glb[:,j]=np.array(da.data[0,:])
            j=j+1
            # print(j)
        else:
            # print(name)
            glb_y[:,jj]=np.array(da.data[-1,:])
            jj=jj+1
            # print(jj)
    else:
        glb_y[:,jj] = np.nansum(np.array(da.data),axis=0)
        jj=jj+1
        # print(jj)

#%% combine with regional:
names_tot=[name for name in names2]
names_tot.extend(names[0:2])
names_y = names[2:len(names)]
glb_mm_month=np.zeros((len(time1),len(names_tot)))
# glb_mm_month[0:len_reg,:]=glb_ts_mm
glb_mm_month[:,0:len_reg]=np.array(df_glb[names2])
glb_mm_month[0:len(time2),len(names2)]=glb[:,0]
glb_mm_month[0:len(time2),len(names2)+1]=glb[:,1]
#%% check time series
for i,n in enumerate(names_tot):
    plt.plot(glb_mm_month[:,i],label=n)
plt.legend()
plt.close()

glb_mm_month[np.where(glb_mm_month==0)]=np.nan
glb_y[np.where(glb_y==0)]=np.nan

# plot again after NaNs
for i,n in enumerate(names_tot):
    plt.plot(glb_mm_month[:,i],label=n)
plt.legend()
plt.close()

#%% create dataset
da = xr.Dataset(data_vars={'OM_monthly_slc':(('months','name_monthly'),glb_mm_month),
                           'OM_yearly_slc':(('years','name_yearly'),-glb_y),
                           
                           },
                         coords={'name_monthly':names_tot,
                                 'name_yearly':names_y,
                                 'months':time1,
                                 'years':years2,
                                 })

da.attrs['units']='mm'
da.attrs['Metadata']='Ocean mean SLC variations from different sources (AIS,GIS,GLA,LWS) and different datasets'
da.attrs['script']='7a-OM_global_source-ts.py'


#%% check dataset
da.OM_monthly_slc.plot(col='name_monthly',col_wrap=4)
plt.close()
da.OM_yearly_slc.plot(col='name_yearly')
plt.close()
#%% save
da.to_netcdf(path+'global_means_timeseries_v2.nc')

path='/Volumes/LaCie_NIOZ/data/barystatic/global/'
da.to_netcdf(path+'global_mean_timeseries_v2.nc')
