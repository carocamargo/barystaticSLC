#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  7 17:31:48 2022

Spatial unc
multiply by mean of period

@author: ccamargo
"""




import sys
sys.path.append("/Users/ccamargo/Documents/py_scripts/")
# import utils_unc as unc
# import utils_SL as sl
import utils_SLE_v2 as sle

import pandas as pd
import xarray as xr
import numpy as np 
# import cmocean as cm
# from cmcrameri import cm as cmf
# import matplotlib.pyplot as plt

#%% open normalized spatial uncertainity:
path = '/Volumes/LaCie_NIOZ/data/barystatic/revisions/results_final/'
file = 'normalized_spatial_unc.p'

df_unc=pd.read_pickle(path+file)
dimlat=180;dimlon=360
X=np.array(df_unc.lon).reshape(dimlat,dimlon)[0,:]
Y=np.array(df_unc.lat).reshape(dimlat,dimlon)[:,0]
#% % open SLF trends

periods =[  # (2005,2016),
          (1993,2018) ,
          (2003,2017)
          ]
for period in periods:
    print(period)
    t0=period[0]
    t1=period[1]-1
    path_to_save = path+'{}-{}/'.format(t0,t1)
    # file = 'SLF_trend_{}-{}.p'.format(t0,t1)
    file = 'SLF_rsl_{}-{}.nc'.format(t0,t1)

    # df_slf = pd.read_pickle(path_to_save+file)
    ds = xr.open_dataset(path_to_save+file)
    
    # one value per dataset
    # df=pd.DataFrame(df_unc[['lon','lat']])
    # for name in df_slf.columns:
    #     Z=np.array(df_slf[name]).reshape(dimlat,dimlon)
    #     _, mu , _ =sle.reg_to_glb(Z,Y,X)
    #     print('{} :{}'.format(name,mu))
    #     reg=name.split('_')[0]
    #     df[name]=df_unc[reg]*mu
        
    # df.to_pickle(path+"spatial_unc_{}-{}.p".format(t0,t1))
    
    # % %open value per region:
    #1.  compute mean trend for each region
    mu = [sle.reg_to_glb(np.array(ds['SLF_trend'][i]),Y,X)[1] for i in range(len(ds.name))]
    # mu = [sle.reg_to_glb(np.array(df_slf[name]).reshape(dimlat,dimlon),Y,X)[1] for name in df_slf.columns]
    reg = [name.split('_')[0] for name in np.array(ds.name)]
    df_mu=pd.DataFrame(mu,columns=['mu'])
    df_mu['reg']=reg
    ens_mean = df_mu.groupby('reg').mean()
    #% %
    # 2. multiply mean by spatial unc:
    df=pd.DataFrame(df_unc[['lon','lat']])
    for region in list(set(reg)):
        print(region)
        df[region] = df_unc[region] * ens_mean.loc[region][0]
    #% %
    df.to_pickle(path_to_save+"spatial_unc_{}-{}.p".format(t0,t1))

#%%
