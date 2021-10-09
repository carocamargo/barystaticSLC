#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 21 13:14:47 2021

Structural ocean mean unc

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
col_dict={1:"black", # WN
          2:"palegoldenrod", # PL
          3:"lightpink", # PLWN
          4:"orange", # AR1
          5:"teal", # Ar5
          6:"darkmagenta", # AR9
          7:"skyblue", # ARf
          8:"crimson" # GGM
          }
def interpolate_gaps(values, limit=None):
    """
    Fill gaps using linear interpolation, optionally only fill gaps up to a
    size of `limit`.
    """
    values = np.asarray(values)
    i = np.arange(values.size)
    valid = np.isfinite(values)
    filled = np.interp(i, i[valid], values[valid])

    if limit is not None:
        invalid = ~valid
        for n in range(1, limit+1):
            invalid[:-n] &= invalid[n:]
        filled[invalid] = np.nan

    return filled
#%%
periods =[ (2005,2016),
            (1993,2018),
            (1993,2017),
            (2003,2017)
          ]
period=periods[0]
for period in periods:
    t0=period[0]
    t1=period[1]
    print('{} to {}'.format(t0,t1-1))

    ifolder = 0
    # t0=2005
    # t1=2016 # until december of previous year
    # pwd='/export/lv1/user/ccamargo/OM/'
    pwd='/Volumes/LaCie_NIOZ/data/barystatic/global/'
    # dataset='means_EWH.nc'
    dataset='global_mean_timeseries_v2.nc'
    #dataset='3-yearly_SLF_monthly_1993-2020_180x360.nc'


    #% % load dataset
    path=pwd
    ds=xr.open_dataset(path+dataset)
    
    #% %
    # ds=da
    # sel time period 
    ds= ds.sel(years=slice(t0,t1-1))
    ds['months'],_=sl.get_dec_time(np.array(ds.months))
    ds= ds.sel(months=slice(t0,t1))
    
    # get variables
    variables = list(ds.keys())
    datasets=[name for name in np.array(ds.name_monthly)]
    if t0<2002:
        datasets = [name for name in np.array(ds.name_monthly) 
                    if not name.endswith('CSR') and not name.endswith('JPL')]
    datasets.extend(name for name in np.array(ds.name_yearly))
    
    dm = pd.DataFrame(np.round(np.array(ds.months),2),columns=['time'])
    dy = pd.DataFrame(np.array(ds.years),columns=['time'])
    
    # ref period:
    tr1 = 2005
    tr2=2015
    ind1_m = np.where(dm['time']==tr1)[0][0]
    ind2_m = np.where(dm['time']==tr2)[0][0]
    ind1_y = np.where(dy['time']==tr1)[0][0]
    ind2_y = np.where(dy['time']==tr2)[0][0]
        
    for iname, name in enumerate(datasets):
        print(name)
        # % %
        if name in np.array(ds.name_yearly):
            idx=np.where(name==np.array(ds.name_yearly))[0][0]
            data = np.array(ds.OM_yearly_slc[:,idx])
            mu = np.nanmean(data[ind1_y:ind2_y])
            dy[name]=np.array(data-mu)


        elif name in np.array(ds.name_monthly):
            idx=np.where(name==np.array(ds.name_monthly))[0][0]
            data = np.array(ds.OM_monthly_slc[:,idx])
            mu = np.nanmean(data[ind1_m:ind2_m])
            dm[name]=data-mu
         
    df = dm.merge(dy,on='time',how='left')
    for name in np.array(dy.columns[1:len(dy.columns)]):
        df[name]=interpolate_gaps(df[name])
        #% %
    regs = ['AIS','GIS','GLA','LWS']
    df_struct=pd.DataFrame(np.array(df['time']),columns=['time'])
    for reg in regs: 
        names = [name for name in datasets if name.startswith(reg)]
        fig=plt.figure()
        df[names].plot()
        #  
        df_struct[reg]=np.array(df[names].std(axis=1))
        df_struct[reg].plot()
        plt.close()
    df_struct[regs].plot()
    df.to_pickle(pwd+'structural_timeseries_{}-{}._v2p'.format(t0,t1-1))
    #% %
    trends = np.zeros((len(datasets)))
    tdec=np.array(df.time)
    for i,name in enumerate(df.columns[1:len(df.columns)]):
        trends[i],error,acc ,_,_= sl.get_OLS_trend(tdec,np.array(df[name]))
        print('{}: {:.3f} ± {:.3f}'.format(name,trends[i],error))
        
    df_tr=pd.DataFrame(trends)
    df_tr['dataset']=np.array(df.columns[1:len(df.columns)])
    df_tr['region'] = [reg.split('_')[0] for reg in df.columns[1:len(df.columns)]]
    df_std = df_tr.groupby(['region']).std()
    df_std = df_std.rename(columns={0:'tr_std'})
    df_std.to_pickle(pwd+'structural_OLS_trend_{}-{}_v2.p'.format(t0,t1-1))
    
    #% % Open hector trends:
    ds=xr.open_dataset(pwd+'hector/ALL-ocean_mean_8_NM-{}-{}_v2.nc'.format(t0,t1-1))
    ic_idx = -1 # BIC_tp
    df_tr=pd.DataFrame(np.array(ds.best_tr[:,ic_idx]))
    df_tr['dataset']=np.array(ds.name)
    df_tr['region'] = [reg.split('_')[0] for reg in df_tr['dataset']]
    df_std = df_tr.groupby(['region']).std()
    df_std = df_std.rename(columns={0:'tr_std'})
    df_std.to_pickle(pwd+'structural_HEC_trend_{}-{}_v2.p'.format(t0,t1-1))
