#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 20 15:08:52 2021

OLS trends

@author: ccamargo
"""



import numpy as np
# import scipy.optimize as opti
import xarray as xr
import matplotlib.pyplot as plt
import sys
sys.path.append("/Users/ccamargo/Documents/py_scripts/")
import utils_SL as sl 
import utils_SLE_v2 as sle

# from netCDF4 import Dataset
import pandas as pd
import os

import datetime as dt

import cmocean as cm
# from mpl_toolkits.basemap import Basemap
# from matplotlib.gridspec import GridSpec
from cartopy import crs as ccrs#, feature as cfeature

#% % packages for plotting
from pandas.plotting import table 
from matplotlib.gridspec import GridSpec
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import ListedColormap

# Let's also design our color mapping: 1s should be plotted in blue, 2s in red, etc...
col_dict={1:"black", # WN
          2:"palegoldenrod", # PL
          3:"lightpink", # PLWN
          4:"orange", # AR1
          5:"teal", # Ar5
          6:"darkmagenta", # AR9
          7:"skyblue", # ARf
          8:"crimson" # GGM
          }

# We create a colormar from our list of colors
cmapnm = ListedColormap([col_dict[x] for x in col_dict.keys()])
def lat2str(deg):
    # Source: https://github.com/matplotlib/basemap/blob/master/examples/customticks.py
    # Adapted so that 0 has no indication of direction.
    minn = 60 * (deg - np.floor(deg)) # transform to minutes
    deg = np.floor(deg) # degrees
    dirr = 'N'
    if deg < 0:
        if minn != 0.0:
            deg += 1.0
            minn -= 60.0
        dirr = 'S'
    elif deg == 0: dirr = ''
    return ("%d\N{DEGREE SIGN} %s") % (np.abs(deg),dirr)

def lon2str(deg):
    # Source: https://github.com/matplotlib/basemap/blob/master/examples/customticks.py
    # Adapted so that 0 has no indication of direction.
    minn = 60 * (deg - np.floor(deg))
    deg = np.floor(deg)
    dirr = ''#'E'
    if deg < 0:
        if minn != 0.0:
            deg += 1.0
            minn -= 60.0
        dirr = ''#'W'
    elif deg == 0: dirr =''
    return ("%d\N{DEGREE SIGN} %s") % (np.abs(deg),dirr)
#% %
# ds_mask = xr.open_dataset('/Volumes/LaCie_NIOZ/data/ETOPO/ETOPO1_Ice-180x360.nc')
# ds_mask.z.plot(# vmin=-1,vmax=1
#                 );#plt.show()

# # ds_mask=ds_mask.sortby('lat',ascending=False)
# oceanmask=np.array(ds_mask.z)
# oceanmask[oceanmask>=0]=1
# oceanmask[oceanmask<=0]=np.nan
#% % 
def ocean_mean(value,lat,lon):
    # value=np.array(ds.best_trend[0,:,:])
    ocean_lit=np.array([360000000,361060000,357000000,360008310,357000000])
    ocean_area= np.mean(ocean_lit)/10**5
    grid_area=sl.get_grid_area(np.ones((180,360)))
    
    
    # plt.pcolor(oceanmask);#plt.show()
    
    # value=np.array(tws_gbw)
    # tdec=np.array(tdec_gwb)
    da=xr.Dataset(data_vars={'data':(('lat','lon'),value)},
                                 coords={
                                         'lat':lat,
                                         'lon':lon})
    mu=(da.data*grid_area).sum(dim=('lat','lon')) /ocean_area
    return mu.data

import pickle
def load_dict(name, path ):
    with open(path + name + '.pkl', 'rb') as f:
        return pickle.load(f)

name = 'masks_dict'
path = '/Volumes/LaCie_NIOZ/data/barystatic/'

masks = load_dict(name,path)
print(masks.keys())

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

  
#%% get all datasets, do OLS trend and unc, run SLE and save 
path_to_save = '/Volumes/LaCie_NIOZ/data/barystatic/OLS/'
path_in='/Volumes/LaCie_NIOZ/data/barystatic/source_var/'
    
dataset_glb='means_v2.nc'
dataset_reg='regional.nc'

# period:
periods = [(2003,2017)]
periods =[ (2005,2016),
          (1993,2018),
          (2003,2017)
          ]
for period in periods:
    t0=period[0]
    t1=period[1]# -1

# period=periods[0]
# t0=period[0]
# t1=period[1]
    print('{} to {}'.format(t0,t1-1))

    #% % 1. OLS on mean:
    ds=xr.open_dataset(path_in+dataset_glb)
    #% % means
    # Select time
    ds= ds.sel(year=slice(t0,t1-1))
    ds= ds.sel(time=slice(t0,t1))

    # get variables
    datasets = list(ds.keys())
    for _, dataset in enumerate(datasets):
        name=dataset
        print(name)

            
        regions = np.array(ds['reg_{}_{}'.format(dataset.split('_')[0],dataset.split('_')[1])])
        tr = np.zeros((len(regions)))
        tr_err = np.full_like(tr, 0)
        acc = np.full_like(tr, 0)


        if dataset.split('_')[1] == 'UCI' or dataset.split('_')[1] == 'ZMP':
            time=np.array(ds.year)
            # sp=365
        else:
            time=np.array(ds.time)
            # sp=30
        data=np.array(ds[dataset])
        # if dataset=='GLA_ZMP_slc':
        for ireg,reg in enumerate(regions):
            height=data[ireg,:]
            acc_idx=np.isfinite(height)
            tr[ireg],tr_err[ireg],acc[ireg] = sl.get_OLS_trend(time[acc_idx], height[acc_idx],ci=90)
     
                
                    
            
        
            da=xr.Dataset(data_vars={'trend':(('reg'),tr),
                                     'unc':(('reg'),tr_err),
                                     'acc':(('reg'),acc),

                                    },
                                    coords={
                                            'fname':dataset,
                                            'reg':regions })
        da.attrs['metadata']='barystatic sea-level trends in mm/y from {} to {}, obtained with OLS model'.format(t0+1,t1)
        da.attrs['uncertainty']='Uncertainty represents the 90% CI of the model'
        da.to_netcdf(path_to_save+'{}_OLS'.format(dataset)+'.nc')

    # % % 2. OLS on regional datasets :
    ds=xr.open_dataset(path_in+dataset_reg)
    #% % means
    # Select time
    ds= ds.sel(year=slice(t0,t1-1))
    to=str(t0)+'-01-01';ti=str(t1)+'-01-01'
    ds['tdec'] = (('time'), ds.tdec)
    ds= ds.sel(time=slice(to,ti))
    # get variables
    datasets = np.array(ds.name)
    lon=np.array(ds.lon)
    lat=np.array(ds.lat)
    time=np.array(ds.tdec)
    i=0;dataset=datasets[0]
    #% % 
    for i, dataset in enumerate(datasets):

        height=np.array(ds.SL_EWH[i,:,:,:])
        tr,tr_err,acc = sl.get_OLS_trend(time, height,lat,lon,ci=90)
        da=xr.Dataset(data_vars={'trend':(('lat','lon'),tr),
                                     'unc':(('lat','lon'),tr_err),
                                     'acc':(('lat','lon'),acc),

                                    },
                                    coords={
                                            'fname':dataset,
                                            'lat':lat,
                                            'lon':lon})
        da.attrs['metadata']='barystatic sea-level trends in mm/y from {} to {}, obtained with OLS model'.format(t0+1,t1)
        da.attrs['uncertainty']='Uncertainty represents the 90% CI of the model'
        da.to_netcdf(path_to_save+'{}_OLS'.format(dataset)+'.nc')

#% % glb to regional
# periods =[
# (2005,2016),
# (1993,2018),
# (2003,2017)
#           ]
# for period in periods:
#     t0=period[0]
#     t1=period[1]-1
    
    # pwd = '/Volumes/LaCie_NIOZ/data/barystatic/hector/source/EWH/'
    pwd = path_to_save
    
    # t0=2005
    # t1=2015
    # pwd = pwd+'{}-{}/'.format(t0,t1)

    flist=sl.get_filelist(pwd,'*.nc')
    # remove intrinsic unc:
    flist = [file for file in flist if not file.split('_')[3]=='unc' and not 
             file.split('/')[-1].split('_')[0]=='TCWS'] 
    # ds=xr.open_dataset(flist[0])
    # nm=np.array(ds.nm)
    # table_tr = np.zeros((len(nm),len(flist)))
    # table_unc = np.full_like(table_tr,0)
    # best_tr = np.zeros((len(flist)))
    # best_unc = np.full_like(best_tr,0)
    name=[None]*len(flist)
    # ic_idx = -1 # use bic_tp
    ifile=1;file=flist[ifile]
    mask_keys = {'AIS_IMB':'AIS_regions',
                 'AIS_UCI':'AIS_basins',
                 'GIS_IMB':'GIS_regions',
                 'GIS_UCI':'GIS_basins',
                 'GLA_ZMP':'Glaciers'
                 }
    dimlat=180;dimlon=360
    for ifile, file in enumerate(flist):
        ds=xr.open_dataset(file)
        name = file.split('/')[-1].split('.')[0]
        if name.split('_')[1]=='300':
            name = name.split('_')[0]+'_'+name.split('_')[2]
        elif name.split('_')[0]=='GLWS':
            name = 'GLA_'+name.split('_')[1]
        
        elif len(name.split('_'))>2:
            name = name.split('_')[0]+'_'+name.split('_')[1]
        
        # name = name.split('_')[0]+'_'+name.split('_')[1]
        print(name)
        if len(ds.trend.shape) ==1: # mean data (reg)
            print(name)
    
            regional_trend = np.zeros((len(ds.reg),dimlat,dimlon))
            regional_unc = np.full_like(regional_trend,0)
            
            ic_idx = 0    
            for ireg,reg in enumerate(np.array(ds.reg)):
                print(reg)
                if name =='GLA_ZMP':
                    reg = reg.split('_')[1]
                mask = masks[mask_keys[name]][reg]
            
                regional_trend[ireg,:,:] = sle.height_to_EWH(glb_to_reg(
                                            np.array(ds.trend[ireg]),mask) ).reshape(180,360) 
                regional_unc[ireg,:,:] = sle.height_to_EWH(glb_to_reg(
                                            np.array(ds.unc[ireg]),mask)  ).reshape(180,360)

                    
            # make regional dataset:
            da=xr.Dataset(data_vars={
                                     'trend':(('reg','lat','lon'),regional_trend),
                                     'unc':(('reg','lat','lon'),regional_unc),
                                     
                       
                            },
                              coords={'lon':np.arange(0.5,360,1),
                                      'lat':np.arange(-89.5,90,1),
                                      'reg':[str(reg) for reg in np.array(ds.reg)],
                                      # 'ic':np.array(ds.ic),
                                      # 'nm':np.array(ds.nm)
                                      })
            if name =='AIS_IMB': 
                da=da.drop_sel(reg='AIS')
                
            elif name =='GIS_IMB' or name =='GIS_UCI':
                da=da.drop_sel(reg='GIS')
            if name.split('_')[0]=='AIS':
                mask=masks['AIS_regions']['AIS']
            elif name.split('_')[0]=='GIS':
                mask=masks['GIS_regions']['GIS']
            else: 
                mask_dic = masks['Glaciers']
                mask_gla=np.zeros((180,360))
                for key in mask_dic:
                    mask_tmp = mask_dic[key]
                    mask_tmp[np.isnan(mask_tmp)]=0
                    mask_gla=mask_gla+mask_tmp
                mask_gla[mask_gla==0]=np.nan
                mask=mask_gla
            da=da.sum(dim='reg')* mask
        else:
            da = ds
        # save 
        
        da.to_netcdf(pwd+'/regional/'+name+'.nc')
            
    #% % open all regional datasets for this time period
    flist=sl.get_filelist(pwd+'regional/','*.nc')
    ds=xr.open_dataset(flist[0])
    trend = np.zeros((len(flist),len(ds.lat),len(ds.lon)))
    unc = np.full_like(trend,0)
    

    name = [file.split('/')[-1].split('.')[0] for file in flist]
    # print(name)
    for i,file in enumerate(flist):
        ds=xr.open_dataset(file)

        trend[i,:,:] = np.array(ds.trend)
        unc[i,:,:] = np.array(ds.unc)
    #% %
    da=xr.Dataset(data_vars={'trend':(('name','lat','lon'),trend),
                             'unc':(('name','lat','lon'),unc),
                                     
                       
                            },
                              coords={'lon':np.array(ds.lon),
                                      'lat':np.array(ds.lat),
                                      'name':name,
                                      # 'ic':np.array(ds.ic),
                                      # 'nm':np.array(ds.nm)
                                      })
    path = '/Volumes/LaCie_NIOZ/data/barystatic/OLS/comb/'
    da.to_netcdf(path+'source_OLS_trend_unc_{}-{}.nc'.format(t0,t1))
        