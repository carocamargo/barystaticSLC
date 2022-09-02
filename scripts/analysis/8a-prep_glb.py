#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 21:31:12 2022

Make gmsl time series

@author: ccamargo
"""


import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
#% %
def make_date(y0,m0,y1,m1):
    date_year=np.arange(y0,y1+1)
    if m1==12:
        y1=y1+1
        m1=1
    else:
        m1 = m1+1
    date_month = np.arange(np.datetime64('{}-{}-15'.format(y0,str(m0).zfill(2))),
                   np.datetime64('{}-{}-15'.format(y1,str(m1).zfill(2))),
        np.timedelta64(1,'M'),dtype='datetime64[M]')
    return date_month , date_year




#%%
def make():
    #% %
    path='/Volumes/LaCie_NIOZ/data/barystatic/revisions/input_data/'
    name = 'regional_v1.nc'
    ds = xr.open_dataset(path+name)
    
    # name = 'means_v1'
    # ds= xr.open_dataset(path+name+'.nc')
    name_save = 'means_reg_v1'
    ds2 = xr.open_dataset(path+name_save+'.nc')
    ds2['time'], _ = make_date(1992,1,2018,12)
    ds2 = ds2.sel(time=slice("1993-01-01", ds.time[-1]))
    
    ds3 = xr.open_dataset(path+'zemp_month_bygla.nc')
    ds3=ds3.sum(dim='reg')
    ds3 = ds3.sel(time=slice("1993-01-01", ds3.time[-1]))
    names = [name for name in np.array(ds.name)]
    names.extend([name for name in np.array(ds2.name)])
    glb = np.zeros((len(names),len(ds.time)))
    glb.fill(np.nan)
    
    for i in range(len(ds.name)):
        glb[i] = np.array(ds['SL_mm'][i].sum(dim=('lat','lon')))
        if ds.name[i]=='LWS_GWB':
            glb[i] = glb[i] * 10 
    for j in range(len(ds2.name)):
        if ds2.name[j]=='GLA_ZMP':
            glb[i+j+1,0:len(ds3.time)] = np.array(ds3['mchange_SLH'])
        else:
            glb[i+j+1,0:len(ds2.time)] = np.array(ds2['SLC_mm'][j].sum(dim=('lat','lon')))
    
    # ds2['SL_EWH'] = ds2['SLC_ewh']
    # da = xr.merge([ds,ds2])
    
    # da = da.sel(time=slice("1993-01-01", ds.time[-1]))
    # ds_glb = da.sum(dim=('lat','lon'))
    #% %
    # y = np.array(ds_glb['SL_mm'])
    gmsl = np.full_like(glb,0)
    for i, name in enumerate(names):

        yy = np.array(glb[i,:])
        # print(y)
        yy[yy==0] = np.nan
        # yy = yy - np.nanmean(yy[120:192]) # 2003-2008
        yy = yy - np.nanmean(yy[120:240]) # 2003-2012
        

        if name.split('_')[1]!= 'IMB':
            gmsl[i,:] = np.array(- yy)
        else:
            gmsl[i,:] = np.array(yy)
    #     plt.plot(gmsl[i],label=name)
    # plt.legend()
    # plt.show()
        #% %
    # intrisic unc time series
    da=xr.open_dataset(path+'intrinsic/SL_unc_JPL_IMBIE.nc')
    names2 = np.array(da.name)
    da = da['UNC_mm'].sum(dim=('lat','lon'))
    
    gmsl_intr = np.array(da.data)
    dx = xr.Dataset(data_vars={'gmsl':(('names','time'),gmsl),
                               'intrisic':(('names2','time'),gmsl_intr),
                               },
                    coords = {'names':names,
                              'names2':names2,
                              'time':ds.time})
    # da.attrs['global_mean'] = 'Global mean removed form 2003-2008'
    dx.attrs['global_mean'] = 'Global mean removed from 2003-2012'
    dx.attrs['units'] = 'mm'
    
    
    dx.to_netcdf(path+'gmsl.nc')
    path = '/Volumes/LaCie_NIOZ/data/barystatic/revisions/results_final/'
    dx.to_netcdf(path+'gmsl.nc')
    
#%%

    
    
#%% run
make()

