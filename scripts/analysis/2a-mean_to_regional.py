#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  7 16:05:06 2022

Redistribute region means into a 1 degree map
@author: ccamargo
"""




import numpy as np
# import scipy.optimize as opti
import xarray as xr
# import matplotlib.pyplot as plt
import sys
sys.path.append("/Users/ccamargo/Documents/py_scripts/")
import utils_SL as sl 
import utils_SLE_v2 as sle

# from netCDF4 import Dataset
import pandas as pd
# import os

# import datetime as dt

# import cmocean as cm
# from mpl_toolkits.basemap import Basemap
# from matplotlib.gridspec import GridSpec
# from cartopy import crs as ccrs#, feature as cfeature

#% % packages for plotting
# from pandas.plotting import table 
# from matplotlib.gridspec import GridSpec
# from mpl_toolkits.basemap import Basemap
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
#                 );##plt.show()

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

def glb_to_reg2(value=1,mask = np.ones((180,360))):
    
    df=pd.DataFrame(mask.flatten(),columns=['mask'])
    df['area'] = sl.get_grid_area(mask).flatten()
    
    df['regional_value'] = (np.full_like(mask,value).flatten() * df['mask'])# df['area']
    # df['regional_value'] = (np.full_like(mask,value).flatten() * df['mask'])/  df['area']

    # df=pd.DataFrame(mask.flatten(),columns=['mask'])
    # df['area'] = sl.get_grid_area(mask).flatten()
    # # df['area']=df['mask']
    # df['weighted_area']=(df['area']*df['mask'])/np.nansum(df['mask']*df['area'])
    # df['regional_value'] = (np.full_like(mask,value).flatten() * df['weighted_area'])
    
    return np.array(df['regional_value']).reshape(mask.shape)

#%% open mask
import pickle
def load_dict(name, path ):
    with open(path + name + '.pkl', 'rb') as f:
        return pickle.load(f)

name_mask = 'masks_dict'
path_mask = '/Volumes/LaCie_NIOZ/data/barystatic/'
# path_mask=path+name
masks = load_dict(name_mask,path_mask)
print(masks.keys())
#%%

periods =[#(2005,2016),
          (1993,2018),
          (2003,2017)
          ]
for period in periods:
    t0=period[0]
    t1=period[1]-1
    
    # pwd = '/Volumes/LaCie_NIOZ/data/barystatic/hector/source/EWH/'
    pwd = '/Volumes/LaCie_NIOZ/data/barystatic/revisions/hector/'
    
    # t0=2005
    # t1=2015
    pwd = pwd+'{}-{}/'.format(t0,t1)

    flist=sl.get_filelist(pwd+'ranking/','*.nc')
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
        if len(ds.trend.shape) ==2: # mean data (reg,nm)
            print(name)
            masks = load_dict(name_mask,path_mask)

            regional_best_trend = np.zeros((len(ds.ic),len(ds.reg),dimlat,dimlon))
            regional_best_trend.fill(np.nan)
            regional_best_unc = np.full_like(regional_best_trend,np.nan)
            regional_rank = np.zeros((len(ds.ic),len(ds.reg),dimlat,dimlon))
            regional_rank.fill(np.nan)
            
            regional_trend = np.zeros((len(ds.reg),len(ds.nm),dimlat,dimlon))
            regional_trend.fill(np.nan)
            regional_unc = np.full_like(regional_trend,np.nan)
            regional_aic = np.full_like(regional_trend,np.nan)
            regional_bic = np.full_like(regional_trend,np.nan)
            regional_bictp = np.full_like(regional_trend,np.nan)
            
            ic_idx = 0    
            for ireg,reg in enumerate(np.array(ds.reg)):
                print(reg)
                if name =='GLA_ZMP':
                    # reg = reg.split('_')[1]
                    regs = ['ALA', 'WNA', 'ACN', 'ACS', 'GRL', 'ISL', 'SJM',
                            'SCA', 'RUA', 'ASN', 'CEU', 'CAU', 'ASC', 'ASW',
                            'ASE', 'TRP', 'SAN', 'NZL', 'ANT']
                    # codes = [i for i in range(len(regs))]
                    # print(regs[int(reg-1)])
                    reg = regs[int(reg-1)]
                    mask = masks[mask_keys[name]][reg]
                    print(reg)
                    for ic_idx in range(len(ds.ic)):
                        # sle.height_to_EWH(data).reshape(180,360)
                        regional_best_trend[ic_idx,ireg,:,:] = glb_to_reg(
                                                np.array(ds.best_trend[ic_idx,ireg]),mask)
                        regional_best_unc[ic_idx,ireg,:,:] = glb_to_reg(
                                                np.array(ds.best_unc[ic_idx,ireg]),mask)  
                        regional_rank[ic_idx,ireg,:,:] = glb_to_reg2(
                                                np.array(ds.ranks[ic_idx,ireg]),mask)   
                        
                    for inm in range(len(ds.nm)):
                        regional_trend[ireg,inm,:,:] = glb_to_reg(
                                                np.array(ds.trend[ireg,inm]),mask) 
                        regional_unc[ireg,inm,:,:] = glb_to_reg(
                                                np.array(ds.unc[ireg,inm]),mask)  
                        regional_aic[ireg,inm,:,:] = glb_to_reg2(
                                                np.array(ds.aic[ireg,inm]),mask)
                        regional_bic[ireg,inm,:,:] = glb_to_reg2(
                                                np.array(ds.aic[ireg,inm]),mask)  
                        regional_bictp[ireg,inm,:,:] =glb_to_reg2(
                                                np.array(ds.aic[ireg,inm]),mask)  
                else:
                        
                    mask = masks[mask_keys[name]][reg]
                    for ic_idx in range(len(ds.ic)):
                        # sle.height_to_EWH(data).reshape(180,360)
                        regional_best_trend[ic_idx,ireg,:,:] = sle.height_to_EWH(glb_to_reg(
                                                np.array(ds.best_trend[ic_idx,ireg]),mask)).reshape(180,360)
                        regional_best_unc[ic_idx,ireg,:,:] = sle.height_to_EWH(glb_to_reg(
                                                np.array(ds.best_unc[ic_idx,ireg]),mask)  ).reshape(180,360)   
                        regional_rank[ic_idx,ireg,:,:] = glb_to_reg2(
                                                np.array(ds.ranks[ic_idx,ireg]),mask)   
                        
                    for inm in range(len(ds.nm)):
                        regional_trend[ireg,inm,:,:] = sle.height_to_EWH(glb_to_reg(
                                                np.array(ds.trend[ireg,inm]),mask) ).reshape(180,360) 
                        regional_unc[ireg,inm,:,:] = sle.height_to_EWH(glb_to_reg(
                                                np.array(ds.unc[ireg,inm]),mask)  ).reshape(180,360)
                        regional_aic[ireg,inm,:,:] = glb_to_reg2(
                                                np.array(ds.aic[ireg,inm]),mask)
                        regional_bic[ireg,inm,:,:] = glb_to_reg2(
                                                np.array(ds.aic[ireg,inm]),mask)  
                        regional_bictp[ireg,inm,:,:] =glb_to_reg2(
                                                np.array(ds.aic[ireg,inm]),mask)  
                        
                        
                    # check : 
                if ic_idx ==0:
                    print('global: {:1.3f}, regional:{:1.3f}'.
                          format(np.array(ds.best_trend[ic_idx,ireg]),
                                 np.nansum(regional_best_trend[ireg,:,:])))
                    
            # make regional dataset:
            da=xr.Dataset(data_vars={'best_trend':(('ic','reg','lat','lon'),regional_best_trend),
                                     'best_unc':(('ic','reg','lat','lon'),regional_best_unc),
                                     'ranks':(('ic','reg','lat','lon'),regional_rank),
                                     
                                     'trend':(('reg','nm','lat','lon'),regional_trend),
                                     'unc':(('reg','nm','lat','lon'),regional_unc),
                                     'aic':(('reg','nm','lat','lon'),regional_aic),
                                     'bic':(('reg','nm','lat','lon'),regional_bic),
                                     'bic_tp':(('reg','nm','lat','lon'),regional_bictp),
                                     
                                     
                       
                            },
                              coords={'lon':np.arange(0.5,360,1),
                                      'lat':np.arange(-89.5,90,1),
                                      'reg':[str(reg) for reg in np.array(ds.reg)],
                                      'ic':np.array(ds.ic),
                                      'nm':np.array(ds.nm)
                                      })
            if name =='AIS_IMB': 
                da=da.drop_sel(reg='AIS')
                
            elif name =='GIS_IMB' or name =='GIS_UCI':
                da=da.drop_sel(reg='GIS')
            if name.split('_')[0]=='AIS':
                mask=masks['AIS_regions']['AIS']
            if name.split('_')[0]=='GIS':
                mask=masks['GIS_regions']['GIS']
            if name=='GLA_ZMP': 
                mask_dic = masks['Glaciers']
                mask_gla=np.zeros((180,360))
                # for key in mask_dic:
                for key in ['ALA', 'WNA', 'ACN', 'ACS', 'ISL', 'SJM',
                        'SCA', 'RUA', 'ASN', 'CEU', 'CAU', 'ASC', 'ASW',
                        'ASE', 'TRP', 'SAN', 'NZL']:
                    print(key)
                    mask_tmp = mask_dic[key]
                    mask_tmp[np.isnan(mask_tmp)]=0
                    mask_gla=mask_gla+mask_tmp
                mask_gla[mask_gla==0]=np.nan
                mask=mask_gla
            da=da.sum(dim='reg')* mask
        else:
            da = ds
        # save 
        
        da.to_netcdf(pwd+'ranking/regional/'+name+'.nc')
            
    #% %
    flist=sl.get_filelist(pwd+'ranking/regional/','*.nc')
    ds=xr.open_dataset(flist[0])
    trend = np.zeros((len(flist),len(ds.nm),len(ds.lat),len(ds.lon)))
    unc = np.full_like(trend,0)
    
    best_trend = np.zeros((len(flist),len(ds.ic),len(ds.lat),len(ds.lon)))
    best_unc = np.full_like(best_trend,0)
    rank = np.full_like(best_trend,0)
    name = [file.split('/')[-1].split('.')[0] for file in flist]
    for i,file in enumerate(flist):
        ds=xr.open_dataset(file)
        best_trend[i,:,:,:] = np.array(ds.best_trend)
        best_unc[i,:,:,:] = np.array(ds.best_unc)
        rank[i,:,:,:] = np.array(ds.ranks)
        trend[i,:,:,:] = np.array(ds.trend)
        unc[i,:,:,:] = np.array(ds.unc)
    #% %
    da=xr.Dataset(data_vars={'best_trend':(('name','ic','lat','lon'),best_trend),
                                     'best_unc':(('name','ic','lat','lon'),best_unc),
                                     'ranks':(('name','ic','lat','lon'),rank),
                                     
                                     'trend':(('name','nm','lat','lon'),trend),
                                     'unc':(('name','nm','lat','lon'),unc),
                                     
                       
                            },
                              coords={'lon':np.array(ds.lon),
                                      'lat':np.array(ds.lat),
                                      'name':name,
                                      'ic':np.array(ds.ic),
                                      'nm':np.array(ds.nm)
                                      })
    path = '/Volumes/LaCie_NIOZ/data/barystatic/revisions/results_final/{}-{}/'.format(t0,t1)
    da.to_netcdf(path+'source_trend_temporal_unc_{}-{}.nc'.format(t0,t1))
