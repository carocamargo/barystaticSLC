#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 29 16:45:35 2021

Run hector for ocean mean

@author: ccamargo
"""


#% % libraries
# hector_reg_bary_9
import numpy as np
import xarray as xr
import sys
sys.path.append("/Users/ccamargo/Documents/py_scripts/")
# import utils_SL as sl 
import utils_hec as hec
import os
# import cmocean as cmo
# from numba import jit


#%%
periods =[ (2005,2016),(1993,2018),
          (2003,2017)]
for period in periods:
    t0=period[0]
    t1=period[1]
    print('{} to {}'.format(t0,t1-1))

    ifolder = 0
    # t0=2005
    # t1=2016 # until december of previous year
    # pwd='/export/lv1/user/ccamargo/OM/'
    pwd='/Volumes/LaCie_NIOZ/data/barystatic/source_var/'
    # dataset='means_EWH.nc'
    dataset='means_v2.nc'
    #dataset='3-yearly_SLF_monthly_1993-2020_180x360.nc'



    #% % set path to hector
    
    path_to_hec='/Volumes/LaCie_NIOZ/dump/'+str(ifolder)+'/'
    os.chdir(path_to_hec)
    #% % load dataset
    path=pwd
    ds=xr.open_dataset(path+dataset)
    
    # path_to_save ='//export/lv1/user/ccamargo/OM/hector/'+str(t0)+'-2016/'
    # path_to_save ='//export/lv1/user/ccamargo/OM/hector/{}-{}/'.format(t0,t1-1)
    # os.system(' cd {}'.format(path_to_save))
    
    # path_to_save ='/Volumes/LaCie_NIOZ/data/barystatic/hector/source/EWH/{}-{}/'.format(t0,t1-1)
    path_to_save ='/Volumes/LaCie_NIOZ/data/barystatic/hector/source/mixed/{}-{}/'.format(t0,t1-1)
    
    #% % means
    factor=10**16
    # Select time
    ds= ds.sel(year=slice(t0,t1-1))
    ds= ds.sel(time=slice(t0,t1))
    
    # get variables
    datasets = list(ds.keys())
    
    # test diff noise models:
    NM=['WN' ,
        'PL',
        'PLWN',
        'AR1',
        'AR5',
        'AR9',
        'ARF',
        'GGMWN',# 'FNWN','RWFNWN'
        ]
    
    
    for _, dataset in enumerate(datasets):
        name=dataset
        print(name)
        if dataset.split('_')[1] == 'UCI' or dataset.split('_')[1] == 'ZMP':
            time=np.array(ds.year)
            sp=365
        else:
            time=np.array(ds.time)
            sp=30
            
        regions = np.array(ds['reg_{}_{}'.format(dataset.split('_')[0],dataset.split('_')[1])])
        tr = np.zeros((len(regions),len(NM)))
        tr_err = np.full_like(tr, 0)
        aic = np.full_like(tr, 0)
        bic = np.full_like(tr, 0)
        bic_c = np.full_like(tr, 0)
        bic_tp = np.full_like(tr, 0)
        logL = np.full_like(tr, 0)
        N = np.full_like(tr, 0)
    
    
        data=np.array(ds[dataset])
        os.chdir(path_to_hec)
        # if dataset=='GLA_ZMP_slc':
        for ireg,reg in enumerate(regions):
            x=data[ireg,:] * factor
            # create a .mom file for it:
            hec.ts_to_mom(x[np.isfinite(x)],time[np.isfinite(x)],
                              sp=sp,path=str(path_to_hec+'raw_files/'),
                              name=str(name),ext='.mom')
            # loop over NM:
            for inm, n in enumerate(NM):
                if n=='ARF':
                    hec.create_estimatetrend_ctl_file(name,n,sp=sp,
                                                          GGM_1mphi = 6.9e-07,
                                                          seas=True,halfseas=True,
                                                          LikelihoodMethod='FullCov'
                                                          )
                else:
                    hec.create_estimatetrend_ctl_file(name,n,sp=sp,
                                                          GGM_1mphi = 6.9e-07,
                                                          seas=False,halfseas=False,
                                                          LikelihoodMethod='FullCov',
                                                          )
    
                # Run estimatetrend (hector)
                os.system('estimatetrend > estimatetrend.out')
    
                # get results:
                out=hec.get_results()
                #save results:
                if out[0]!=None:
                    tr[ireg,inm]=float(out[0])/factor
                    tr_err[ireg,inm]=float(out[1])/factor
                    N[ireg,inm]=out[2]
                    logL[ireg,inm]=out[3]
                    aic[ireg,inm]=out[4]
                    bic[ireg,inm]=out[5]
                    bic_tp[ireg,inm]=out[6]
                    bic_c[ireg,inm]=out[7]
                
                    
            
        
            da=xr.Dataset(data_vars={'trend':(('reg','nm'),tr),
                                     'unc':(('reg','nm'),tr_err),
                                     'aic':(('reg','nm'),aic),
                                     'bic':(('reg','nm'),bic),
                                     'bic_c':(('reg','nm'),bic_c),
                                     'bic_tp':(('reg','nm'),bic_tp),
                                     'logL':(('reg','nm'),logL),
                                     'N':(('reg','nm'),N),
                                    },
                                    coords={'nm':NM,
                                            'fname':dataset,
                                            'reg':regions })
        da.attrs['metadata']='barystatic sea-level trends in mm/y from {} to {}, obtained with Hector'.format(t0+1,t1)
        os.system(' cd {}'.format(path_to_save))
        da.to_netcdf(path_to_save+'{}_{}_NM'.format(dataset,len(NM))+'.nc')

