#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  6 15:06:22 2021

@author: ccamargo
"""





import pandas as pd
import numpy as np 
import xarray as xr
from matplotlib import pyplot as plt
from cartopy import crs as ccrs# , feature as cfeature

def quadrsum(X):
    Y = np.zeros((X.shape[1],X.shape[2]))
    for i in range(X.shape[0]):
        Y = Y + (X[i]*X[i])
    Y = np.sqrt(Y)
    return Y



X=np.random.rand(3,4,5)
Y= quadrsum(X)

Y2 = np.sqrt((X[0]*X[0])+(X[1]*X[1])+(X[2]*X[2]))



# t0=2005
# t1=2015

periods =[
        (2005,2016),
          (1993,2018),
          (2003,2017)
          ]
#%%
for period in periods:
    #% %
    # period=periods[-1]
    t0=period[0]
    t1=period[1]-1
    path='/Volumes/LaCie_NIOZ/data/barystatic/results/{}-{}/'.format(t0,t1)
    df_temporal_unc = pd.read_pickle(path+'SLF_temporal_unc_{}-{}.p'.format(t0,t1))
    
    # ds_trend =xr.open_dataset(path+'SLF_trend_ALL_NM_OLS_{}-{}.nc'.format(t0,t1))
    # ds_trend=ds_trend.sel(nm='OLS')
    # trend=np.array(ds_trend.trend).reshape(len(ds_trend.name),len(ds_trend.lat)*len(ds_trend.lon))
    # df_trend = pd.DataFrame(trend.T,columns=np.array(ds_trend.name))
    
    df_trend = pd.read_pickle(path+'SLF_trend_{}-{}.p'.format(t0,t1))
    
    # ds_intrinsic_unc = xr.open_dataset(path+'intrinsic_unc_{}-{}.nc'.format(t0,t1))
    # ds_intrinsic_unc = xr.open_dataset(path+'intrinsic_unc_source_trend_SLF{}-{}.nc'.format(t0,t1))
    ds_intrinsic_unc = xr.open_dataset(path+'intrinsic_unc_prop_{}-{}.nc'.format(t0,t1))
    names = [name for name in np.array(ds_intrinsic_unc.name)]
    for iname, name in enumerate(names):
        if name=='GLWS_JPL': names[iname] = 'GLA_JPL'
    ds_intrinsic_unc['name']=names    
    ds_noise_model = xr.open_dataset(path+'source_trend_temporal_unc_{}-{}.nc'.format(t0,t1))
    # ds_noise_model = ds_noise_model.sortby('lat',ascending=False)
    df_spatial_unc = pd.read_pickle(path+'spatial_unc_{}-{}.p'.format(t0,t1))
    
    path2='/Volumes/LaCie_NIOZ/data/barystatic/results/'
    file='normalized_spatial_unc.p'
    df_spatial_unc_norm=pd.read_pickle(path2+file)
    
    names = [name for name in df_trend.columns]
    lat=np.array(ds_intrinsic_unc.lat)
    lon=np.array(ds_intrinsic_unc.lon)
    dimlat=len(lat)
    dimlon=len(lon)
    trends=np.zeros((len(names),len(lat),len(lon)))
    unc_type=['temporal','spatial','intrinsic']
    uncs = np.zeros((len(unc_type),len(names),len(lat),len(lon)))
    nm_sel = np.full_like(trends,0)
    spatial_unc = np.full_like(trends,0)
    ic_idx = -1
    # this should be the same IC index as in the script 2.e.SLE_OM_hectorsource.py, 
    # when we inputed the best trend and uncertainty in the SLE! 
    
    for iname, name in enumerate(names):
        trends[iname,:,:]=np.array(df_trend[name]).reshape(dimlat,dimlon)
        # trends[iname,:,:]=np.array(ds_noise_model.sel(name=name).best_trend[ic_idx,:,:])
        nm_sel[iname,:,:]=np.array(ds_noise_model.sel(name=name).ranks[ic_idx,:,:])
        
        reg=name.split('_')[0]
        spatial_unc[iname,:,:]=np.array(df_spatial_unc_norm[reg]).reshape(dimlat,dimlon)
        for iunc, unc in enumerate(unc_type):
            if unc=='temporal':
                uncs[iunc,iname,:,:]=np.abs(np.array(df_temporal_unc[name]).reshape(dimlat,dimlon))
                # uncs[iunc,iname,:,:]=np.abs(ds_noise_model.sel(name=name).best_unc[ic_idx,:,:])
                
            elif unc =='spatial':
                uncs[iunc,iname,:,:]=np.abs(np.array(df_spatial_unc[name.split('_')[0]]).reshape(dimlat,dimlon))
            else: # unc=='intrinsic'
                if name in np.array(ds_intrinsic_unc.name):
                    # uncs[iunc,iname,:,:]=np.abs(np.array(ds_intrinsic_unc.sel(name=name).intrinsic_unc_SLF_trend))
                    uncs[iunc,iname,:,:]=np.abs(np.array(ds_intrinsic_unc.sel(name=name).intrinsic_unc_SLF))
                
    #% % sum the uncertainities in quadrature:
    unc_total = np.sqrt((uncs[0]**2)+
                        (uncs[2]**2)+
                        (uncs[1]**2))
    # unc_total = uncs.sum(axis=0)        
    
    #% % make dataset
    da=xr.Dataset(data_vars={'trend':(('name','lat','lon'),trends),
                              'uncs':(('unc_type','name','lat','lon'),uncs),
                              'unc_total':(('name','lat','lon'),unc_total),
                              'nm_sel':(('name','lat','lon'),nm_sel),
                              'spatial_unc_norm':(('name','lat','lon'),spatial_unc),
                             
    
                            },                          coords={'lon':lon,
                                      'lat':lat,
                                      'name':names,
                                      'unc_type':unc_type,
                                      'nm':np.array(ds_noise_model.nm)
                                      })
                                                                
    # da.uncs[:,0,:,:].plot(col='unc_type',cmap='RdBu_r',vmin=-1,vmax=1)
    # plt.show()
    #% % save dataset
    da.to_netcdf(path+'final_dataset_OLS-prop_{}-{}.nc'.format(t0,t1))
    
    # % make and save combos
    # ds=xr.open_dataset(path+'final_dataset_{}-{}.nc'.format(t0,t1))
    ds= da
    ds=ds.sortby('lat',ascending=False)
    
    lon=np.array(ds.lon)
    lat=np.array(ds.lat)
    llon,llat=np.meshgrid(lon,-lat)
    df=pd.DataFrame(np.hstack(llon),columns=['lon'])
    df['lat']=np.hstack(llat)
    
    #% %
    combos = [['JPL'],
                ['CSR'],
               ['IMB','WGP'],
              ['IMB','GWB','ZMP'],['UCI','WGP'],['UCI','GWB','ZMP']
             ]
    reconstr =['JPL',
                'CSR',
                'IMB+WGP',
               'IMB+GWB+ZMP','UCI+WGP','UCI+GWB+ZMP'
              ]
    # i=0
    for title,combo in zip(reconstr,combos):
        
        #% %
        names=[name for name in np.array(ds.name) for comb in combo if name.split('_')[1]==comb ]
        regions = [name.split('_')[0] for name in np.array(ds.name) for comb in combo if name.split('_')[1]==comb]
        # trend=np.nansum(np.array(ds.trend.sel(name=names)),axis=0)    
        # unc=quadrsum(np.array(ds.unc_total.sel(name=names)))
        df['{}_trend_tot'.format(title)]=np.hstack(np.nansum(np.array(ds.trend.sel(name=names)),axis=0))
        df['{}_unc_tot'.format(title)]=np.hstack(quadrsum(np.array(ds.unc_total.sel(name=names))))
        # da=ds.sel(name=names)
        for i,unc_typ in enumerate(np.array(ds.unc_type)):
            # print(unc_typ)
            df['{}_unc_{}'.format(title,unc_typ)]=quadrsum(np.array(ds.uncs[i].sel(name=names))).flatten()
        for i,reg in enumerate(regions):
            df['{}_unc_{}'.format(title,reg)]=quadrsum(np.array(ds.uncs.sel(name=names[i]))).flatten()
        
        df.to_pickle(path+'OM_reconstructions_OLS-prop_{}-{}.p'.format(t0,t1))
        
        
        
        