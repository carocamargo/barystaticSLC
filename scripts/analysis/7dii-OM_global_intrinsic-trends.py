#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 22 13:13:38 2021

Global intrisinc unc trend

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
#%% open dataset
path='/Volumes/LaCie_NIOZ/data/barystatic/global/'

periods =[ (2005,2016),
            (1993,2018),
            (1993,2017),
            (2003,2017)
          ]

# loop over time: 
period=periods[0]
for period in periods:
    ds=xr.open_dataset(path+'global_intrinsic-timeseries_v2.nc')
    ds_sl=xr.open_dataset(path+'global_mean_timeseries_v2.nc')
    
    t0=period[0]
    t1=period[1]
    print('{} to {}'.format(t0,t1-1))


    datasets=[name for name in np.array(ds.name_monthly)]
    if t0<2002:
        datasets = [name for name in np.array(ds.name_monthly) 
                    if not name.endswith('CSR') and not name.endswith('JPL')]
    

    # select time period
    ds['months'],_=sl.get_dec_time(np.array(ds.months))
    ds= ds.sel(months=slice(t0,t1))
    ds_sl['months'],_=sl.get_dec_time(np.array(ds_sl.months))
    ds_sl= ds_sl.sel(months=slice(t0,t1))
    
    trends = np.zeros((len(datasets)))
    std_trend = np.full_like(trends, 0 )
    tdec= np.array(ds.months)
    for i,name in enumerate(datasets):
        da=ds.sel(name_monthly=name)
        da_sl=ds_sl.sel(name_monthly=name)
        
        sigma = np.array(da.OM_monthly_intrinsic_unc )
        y = np.array(da_sl.OM_monthly_slc )
        
        trends[i],error,acc, trend2, std_trend[i] = sl.get_OLS_trend(tdec,y, sigma=sigma)
        print('{}: {:.3f} ± {:.3f}'.format(name,trends[i],std_trend[i]))
    
    df_tr=pd.DataFrame(np.abs(std_trend)) # unc needs to be positive
    df_tr['dataset']=datasets
    df_tr = df_tr.rename(columns={0:'trend'})
    df_tr.to_pickle(path+'intrinsic_OLS_trend_{}-{}_v2.p'.format(t0,t1-1))
    
    #%%
    