#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  2 17:24:58 2021

@author: ccamargo
"""

#%% libraries
# hector_reg_bary_9
import numpy as np
import xarray as xr
import sys
sys.path.append("/export/lv1/user/ccamargo/py_scripts/")
# import utils_SL as sl 
import utils_hec as hec
import os
# import cmocean as cmo
# from numba import jit
import datetime as dt

#%% open dataset
pwd='/export/lv1/user/ccamargo/OM/'
dataset='12_input_datasets_buf_1993-2020_180x360_update_v2.nc'
ds = xr.open_dataset(pwd+dataset)
datasets = np.array(ds.name)
#set decimal time as variable
ds['tyear']=(('time'),ds.tdec)
lat=np.array(ds.lat)
lon=np.array(ds.lon)
dimlat=len(lat);dimlon=len(lon)
lon,lat=np.meshgrid(lon,lat)
lon=np.hstack(lon)
# lon3=lon.flatten() # same thing
lat=np.hstack(lat) # or lat.flatten()

# test diff noise models:
NM=['WN','PL','PLWN','AR1','AR5','AR9','ARF','GGMWN',# 'FNWN','RWFNWN'
    ]

#%% loop over time
periods =[ 
        (2005,2016),
        (1993,2018),
        (2003,2017)
          ]
for period in periods:
    t0=period[0]
    t1=period[1]
    print('{} to {}'.format(t0,t1-1))

    ifolder = 0    
    path_to_hec='/export/lv1/user/ccamargo/dump/'+str(ifolder)+'/'
    os.chdir(path_to_hec)
    # Select time
    to=str(t0)+'-01-01';ti=str(t1)+'-01-01'
    da= ds.sel(time=slice(to,ti))

    # loop over datasets:
    for iname, name in enumerate(datasets):
        print(str(name))
        # Select dataset
        da=da.sel(name=[name])
        data=np.array(da.SL_mm[iname,:,:,:])
        time=np.array(da.tyear)
        data2=data.ravel().reshape(len(time),len(lon))
        
        # allocate empty variables
        bic=np.full_like(np.zeros((len(NM),len(lon))),np.nan)
        bic_c=np.full_like(np.zeros((len(NM),len(lon))),np.nan)
        bic_tp=np.full_like(np.zeros((len(NM),len(lon))),np.nan)
        aic=np.full_like(np.zeros((len(NM),len(lon))),np.nan)
        logL=np.full_like(np.zeros((len(NM),len(lon))),np.nan)
        N=np.full_like(np.zeros((len(NM),len(lon))),np.nan)
        tr=np.full_like(np.zeros((len(NM),len(lon))),np.nan)
        tr_err=np.full_like(np.zeros((len(NM),len(lon))),np.nan)



        path_to_save ='/export/lv1/user/ccamargo/OM/hector/{}-{}/'.format(t0,t1-1)
        os.system(' cd {}'.format(path_to_save))
        inm=0;n=NM[inm]

        #% % Loop over each lat-lon then over NM
        go = dt.datetime.now()
        for ilon in range(len(lon)):
            if np.any(np.isfinite(data2[:,ilon])): # check if we have data:
                print(ilon)
                x=data2[:,ilon]
                # create a .mom file for it:
                hec.ts_to_mom(x[np.isfinite(x)],time[np.isfinite(x)],
                              sp=31,path=str(path_to_hec+'raw_files/'),
                              name=str(name),ext='.mom')
                
                # loop over NM:
                for inm, n in enumerate(NM):
                    # Create a .ctl file:
                    hec.create_estimatetrend_ctl_file(name,n,sp=31,GGM_1mphi = 6.9e-07,LikelihoodMethod='FullCov')
            
                    # Run estimatetrend (hector)
                    os.system('estimatetrend > estimatetrend.out')
            
                    # get results:
                    out=hec.get_results()
                    #save results:
                    if out[0]!=None:
                        tr[inm,ilon]=out[0]
                        tr_err[inm,ilon]=out[1]
                        N[inm,ilon]=out[2]
                        logL[inm,ilon]=out[3]
                        aic[inm,ilon]=out[4]
                        bic[inm,ilon]=out[5]
                        bic_c[inm,ilon]=out[6]
                        bic_tp[inm,ilon]=out[7]
        
        print( dt.datetime.now() - go)
        # save all
        print('saving for: '+str(name))
        dsh=xr.Dataset(data_vars={'trend':(('nm','lat','lon'),tr.reshape(len(NM),dimlat,dimlon)),
                              'unc':(('nm','lat','lon'),tr_err.reshape(len(NM),dimlat,dimlon)),
                              'bic':(('nm','lat','lon'),bic.reshape(len(NM),dimlat,dimlon)),
                              'aic':(('nm','lat','lon'),aic.reshape(len(NM),dimlat,dimlon)),
                              'logL':(('nm','lat','lon'),logL.reshape(len(NM),dimlat,dimlon)),
                               'bic_c':(('nm','lat','lon'),bic_c.reshape(len(NM),dimlat,dimlon)),
                              'bic_tp':(('nm','lat','lon'),bic_tp.reshape(len(NM),dimlat,dimlon)),
                              'N':(('nm','lat','lon'),N.reshape(len(NM),dimlat,dimlon)),
                               },
                                coords={'nm':NM,
                                        'fname':str(name),                                
                                        'lat':lat.reshape(dimlat,dimlon)[:,0],
                                        'lon':lon.reshape(dimlat,dimlon)[0,:]})
        dsh.attrs['metadata']='barystatic sea-level trends in mm/y from {} to {}, obtained with Hector'.format(
            t0,t1-1)
        dsh.to_netcdf(path_to_save+str(name)+'_{}_NM'.format(len(NM))+'.nc')
        print('saved for: '+str(name))
