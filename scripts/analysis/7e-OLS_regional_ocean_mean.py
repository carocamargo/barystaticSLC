#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 24 09:15:39 2021

Make OLS trends on regional scale, run SLE, 
and then compute ocean mean
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
#%%
import pickle
def load_dict(name, path ):
    with open(path + name + '.pkl', 'rb') as f:
        return pickle.load(f)
def ocean_mean(value,lat,lon):
    # value=np.array(ds.best_trend[0,:,:])
    ocean_lit=np.array([360000000,361060000,357000000,360008310,357000000])
    ocean_area= np.mean(ocean_lit)/10**5
    grid_area=sl.get_grid_area(np.ones((180,360)))
    
    
    # plt.pcolor(oceanmask);##plt.show()
    
    # value=np.array(tws_gbw)
    # tdec=np.array(tdec_gwb)
    da=xr.Dataset(data_vars={'data':(('lat','lon'),value)},
                                 coords={
                                         'lat':lat,
                                         'lon':lon})
    mu=(da.data*grid_area).sum(dim=('lat','lon')) /ocean_area
    return mu.data
def glb_to_reg(value=1,mask = np.ones((180,360))):
    
    df=pd.DataFrame(mask.flatten(),columns=['mask'])
    df['area'] = sl.get_grid_area(mask).flatten()
    
    df['regional_value'] = (np.full_like(mask,value).flatten() * df['mask'])/ len(mask[np.isfinite(mask)])# df['area']
    # df['regional_value'] = (np.full_like(mask,value).flatten() * df['mask'])/  df['area']

    # df=pd.DataFrame(mask.flatten(),columns=['mask'])
    # df['area'] = sl.get_grid_area(mask).flatten()
    # # df['area']=df['mask']
    # df['weighted_area']=(df['area']*df['mask'])/np.nansum(df['mask']*df['area'])
    # df['regional_value'] = (np.full_like(mask,value).flatten() * df['weighted_area'])
    
    return np.array(df['regional_value']).reshape(mask.shape)

#%% loop over time

periods =[ 
            # (2005,2016),
            # (1993,2018),
            (1993,2017),
            (2003,2017)
          ]
period=periods[0]

for period in periods:
    t0=period[0]
    t1=period[1]
     #% %  regional 
    path ='/Volumes/LaCie_NIOZ/data/barystatic/source_var/'
    file1='regional_v3.nc'
    
    ds=xr.open_dataset(path+file1)
    ds=ds.drop_sel(name=['TCWS_CSR','TCWS_JPL'])
    time1=np.array(ds.time)
    tdec,tdec0=sl.get_dec_time(time1)
    ds['time']=tdec
    ds= ds.sel(time=slice(t0,t1))
    tdec=np.array(ds.time)
    
    lat=np.array(ds.lat)
    lon=np.array(ds.lon)
    dimlon=len(lon)
    dimlat=len(lat)
    
    #% % namelist
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
    
    
    #% % compute ocean mean
    iname=0;name=names2[iname]
    df_tr = pd.DataFrame(names2,columns=['dataset'])
    mu=np.zeros((len(names2)))
    df_glb=pd.DataFrame(np.hstack(np.array(tdec)),columns=['time'])
    iname=0;name=names2[iname]
    for iname, name in enumerate(names2):
        ds3 = ds.sel(name=name)
        
        value=np.array(ds3.SL_mm[:,:,:])
        # get regional trend:
        trend, _ , _,_,_ =sl.get_OLS_trend(tdec,value,lat=lat,lon=lon)
        # run SLE:
        slf=sle.run_SLE(sle.height_to_vol(trend).reshape(dimlat,dimlon),name+'_OLS') 
        # compute ocean mean:
        mu[iname],_,_ = sle.reg_to_glb(slf, lat, lon)
        
        
        # out=sle.run_SLE(sle.height_to_vol(out[0]).reshape(dimlat,dimlon),'_test')
    df_tr['trend']=mu
    
    
    #% % basin means
    file2 = 'means_v2.nc'
    ds=xr.open_dataset(path+file2)
    
    ds= ds.sel(time=slice(t0,t1))
    ds= ds.sel(year=slice(t0,t1))
    
    datasets = [dataset for dataset in list(ds.keys()) if dataset.endswith('slc')]
    # reg_dic = {}
    
    masks = load_dict('masks_dict','/Volumes/LaCie_NIOZ/data/barystatic/')
    mask_keys = {'AIS_IMB':'AIS_regions',
                     'AIS_UCI':'AIS_basins',
                     'GIS_IMB':'GIS_regions',
                     'GIS_UCI':'GIS_basins',
                     'GLA_ZMP':'Glaciers'
                     }
    for _, dataset in enumerate(datasets):
        name=dataset.split('_')[0]+'_'+dataset.split('_')[1]
        print(name)
        if dataset.split('_')[1] == 'UCI' or dataset.split('_')[1] == 'ZMP':
            time=np.array(ds.year)
            # sp=365
        else:
            time=np.array(ds.time)
            sp=30
            
        regions = np.array(ds['reg_{}_{}'.format(dataset.split('_')[0],dataset.split('_')[1])])
        tr = np.zeros((len(regions)))
    
        data=np.array(ds[dataset])
    
        # if dataset=='GLA_ZMP_slc':
        for ireg,reg in enumerate(regions):
            x=data[ireg,:] 
            if ~np.all(x==0):
                tr[ireg], _, _,_,_ = sl.get_OLS_trend(time[np.isfinite(x)],x[np.isfinite(x)])
            print(' {} : {:.3f}'.format(reg,tr[ireg]))
        # end for loop of regions
        # reg_dic[name+'_trend']=tr
        # reg_dic[name+'_regions']=regions
        
        # transform from basin mean to regional:
        regional_trend = np.zeros((len(regions),dimlat,dimlon))
        for ireg,reg in enumerate(regions):
            print(reg)
            if name =='GLA_ZMP':
                reg = reg.split('_')[1]
            # get mask to regrid to regional
            mask = masks[mask_keys[name]][reg]
            regional_trend[ireg,:,:] = glb_to_reg(tr[ireg],mask).reshape(dimlat,dimlon)
        # reg_dic[name+'trend_by_reg']=regional_trend
        # combine different regions into a single matrix:
        # remove total AIS, and GIS (to not double count)
        regional_trend[np.isnan(regional_trend)]=0
        if name.split('_')[0] in regions: 
            idx=np.where(name.split('_')[0]!=regions)[0]
            r=regional_trend[idx,:,:].sum(axis=0)
        else:
            r=regional_trend.sum(axis=0)
        r[np.where(r==0)]=np.nan
        # plt.pcolor(r);plt.title(name);plt.show()
        
        # run SLE with the regional trend
        slf=sle.run_SLE(sle.height_to_vol(r).reshape(dimlat,dimlon),name+'_OLS') 
        # compute ocean mean:
        mu,_,_ = sle.reg_to_glb(slf, lat, lon)
        df_tr = df_tr.append({'dataset' : name,
                    'trend' : np.abs(np.array(mu))} , 
                    ignore_index=True)
    #% %
    print(df_tr)
    path='/Volumes/LaCie_NIOZ/data/barystatic/global/'
    df_tr.to_pickle(path+'OLS_SLF_reg_mean-{}-{}.p'.format(t0,t1-1))
