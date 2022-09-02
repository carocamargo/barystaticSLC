#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 18 11:09:27 2022

@author: ccamargo
"""
import numpy as np
import xarray as xr
import sys
sys.path.append("/Users/ccamargo/Documents/py_scripts/")
import utils_SL as sl 
# import utils_SLE_v2 as sle
import utils_hec as hec
sys.path.append("/Users/ccamargo/Documents/py_scripts/OM/")
import utils_OM as om 
# from datetime import datetime as dt
import os

import pandas as pd
#%%
def open_gmsl():
    path = '/Volumes/LaCie_NIOZ/data/barystatic/revisions/results_final/'
    ds = xr.open_dataset(path+'gmsl.nc')
    return ds

def intr_unc(x,y,s):
    return om.intrinsic(x,y,s,1000,dim='1D',ci_level=0.95)

#% %
def find_best_val(tr,unc,idx):
    """
    tr,unc and idx are 1D arrays with the length of the used noise models.
    tr are the trends for each NM, unc the unc for each NM, 
    and idx is the AIC, BIC or BICtp ranking for each noise model
    """
    target = idx # all NM for this time series
    len_nm= len(idx)
    logidx= np.array([np.exp((np.nanmin(target)-target[inm])/2) for inm in range(0,len_nm)])
    
    threshold = 0.5
    if np.any(logidx>threshold):
        ind=np.where(logidx>0.5)[0]
        if len(ind)>1:
            ind=np.where(tr==np.min(tr[ind]))[0][0]
        else:
            ind=ind[0]
    else:
        ind=np.where(np.nanmin(target))[0][0]
    
    best_tr=tr[ind]
    best_unc=unc[ind]
    
    return best_tr,best_unc,ind
#% %
def hector(x,y,nm,sp):
    ''' 
    x = time
    y = variable
    nm = noise model
    sp = sampling frequency
    
    '''

    path_to_hec='/Volumes/LaCie_NIOZ/dump/0/'
    os.chdir(path_to_hec)
    name='gmsl'
    # create a .mom file for it:
    hec.ts_to_mom(y[np.isfinite(y)],x[np.isfinite(y)],
                      sp=sp,path=str(path_to_hec+'raw_files/'),
                      name=str(name),ext='.mom')
    
    hec.create_estimatetrend_ctl_file(name,nm,sp=sp,
                                          GGM_1mphi = 6.9e-07,
                                              seas=True,halfseas=True,
                                              LikelihoodMethod='FullCov',
                                              )
    os.system('estimatetrend > estimatetrend.out')

    # get results:
    out=hec.get_results()
    # trend, unc, N, logL, aic, bic, bic_tp, bic_c = out
    return out[0:8]

def QC(trend, unc, aic, bic, bic_tp):
    std = np.nanstd(trend) + np.nanmax(unc)
    mu = np.nanmean(trend)
    if np.any(trend> mu +std) or np.any(trend< mu-std):
        idx=np.array(np.where((trend > mu+std) | (trend < mu-std)))[0]
        for ix in idx:
            # print('QC removed {}'.format(nm[ix]))
            trend[ix]=np.nan;
            unc[ix]=np.nan
            # aic[jname,ix,ilat,ilon]=np.nan;
            # bic[jname,ix,ilat,ilon]=np.nan;
            # logL[jname,ix,ilat,ilon]=np.nan
            aic[ix]=np.nanmax(aic)
            bic[ix]=np.nanmax(bic)
            bic_tp[ix]=np.nanmax(bic_tp)
            # logL[x]=np.nanmax(logL)
    return [trend, unc, aic, bic, bic_tp]

def hector_loop(x,y, NM=['WN' ,
     'PL',
     'PLWN',
     'AR1',
     'AR5',
     'AR9',
     'ARF',
     'GGMWN',# 'FNWN','RWFNWN'
     # 'OLS'
     ],sp=30):
    
    trends = np.zeros((len(NM)))
    uncs = np.zeros((len(NM)))
    aic = np.zeros((len(NM)))
    bic = np.zeros((len(NM)))
    bic_tp = np.zeros((len(NM)))
    
    for i,n in enumerate(NM):
        trends[i], uncs[i],_, _, aic[i], bic[i], bic_tp[i], _ = hector(x,y,n,sp)
    
    trends, uncs, aic, bic, bic_tp = QC(trends, uncs, aic, bic, bic_tp)
    best_tr = np.zeros((3))
    best_unc = np.zeros((3))
    best_nm = np.zeros((3))
    
    idx = 0
    best_tr[idx],best_unc[idx],best_nm[idx] = find_best_val(trends, uncs, aic)
    idx = 1
    best_tr[idx],best_unc[idx],best_nm[idx] = find_best_val(trends, uncs, bic)
    idx = 2
    best_tr[idx],best_unc[idx],best_nm[idx] = find_best_val(trends, uncs, bic_tp)
    best_NM = [NM[int(i)] for i in best_nm]

    return best_tr[-1],best_unc[-1], best_nm[-1], best_NM[-1]

def structural(trends,names):
    regs = [name.split('_')[0] for name in names]
    df = pd.DataFrame({'trends':trends,
                       'names':names,
                       'region':regs})
    d = df.groupby('region').std()
    df['structural'] = [d.loc[reg][0] for reg in regs]
    return df
#%%
def make():
    
    periods = [
        (2003,2016),
        (1993,2016)
        ]
    for period in periods:
        t0,t1 = period
        print(period)
        #% %
        # t0,t1=(2003,2016)
        ds = open_gmsl()
        ds = ds.sel(time=slice("{}-01-01".format(t0), "{}-12-01".format(t1)))
        x,_ = sl.get_dec_time(np.array(ds.time))
        
        names = np.array(ds.names)
        if t0<2002:
            names = [name for name in names if name.split('_')[1]!='JPL' and name.split("_")[1]!='CSR']
        intr = np.zeros((len(names)))
        temp = np.zeros((len(names)))
        trend = np.zeros((len(names)))
        nm = []
        nm_idx = np.zeros((len(names)))
        
        for iname, name in enumerate(names):
            print(name)
            y = np.array(ds['gmsl'].sel(names=name))
            idx = np.isfinite(y)
            trend[iname],temp[iname],nm_idx[iname], n = hector_loop(x, y)
            nm.append(n)
            if name in np.array(ds.names2):
                s = np.array(ds['intrisic'].sel(names2=name))    
                
                intr[iname] = intr_unc(x[idx], y[idx], s[idx])
        df = structural(trend, names)
        df['temporal']= temp
        df['intrinsic'] = intr
        struc = np.array(df['structural'])
        total_unc = np.array( np.sqrt( (temp**2 + intr**2 + struc**2) )  )
        df['uncertainty'] = total_unc
        
        
        ds['trends'] = (('name'),trend)
        ds['temporal_unc'] = (('name'),temp)
        ds['NM'] = (('name'),nm)
        ds['NM_idx'] = (('name'),nm_idx)
        ds['intrinsic_unc'] = (('name'),intr)
        ds['total_unc'] = (('name'),total_unc)
        
        if t0==1993:
            t1=2017
        path='/Volumes/LaCie_NIOZ/data/barystatic/revisions/results_final/{}-{}/'.format(t0,t1)
        ds.to_netcdf(path+'gmsl_ts_trends.nc')
        df.to_pickle(path+'gmsl_trends.p')
    

#%% run
make()
