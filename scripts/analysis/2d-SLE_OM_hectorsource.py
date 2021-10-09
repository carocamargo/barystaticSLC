#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  2 11:58:27 2021

run SLE with the best trend and unc 

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
#%%
periods =[ (2005,2016),(1993,2018),
          (2003,2017)
          ]
for period in periods:
    t0=period[0]
    t1=period[1]-1
    
    # pwd = '/Volumes/LaCie_NIOZ/data/barystatic/hector/source/EWH/'
    pwd = '/Volumes/LaCie_NIOZ/data/barystatic/hector/source/mixed/'
    # t0=2005
    # t1=2015
    pwd = pwd+'{}-{}/'.format(t0,t1)

    flist=sl.get_filelist(pwd+'ranking/regional/','*.nc')
    
    ic_idx=-1
    for ifile, file in enumerate(flist):
        ds=xr.open_dataset(file)
        name = file.split('/')[-1].split('.')[0]
        dataset = name.split('_')[1]
        # SLE
        X=np.array(ds.lon); Y=np.array(-ds.lat)
        data = np.array(ds.best_trend[ic_idx,:,:])
        if dataset=='IMB':
            data=-data
        # if name[ifile].split('_')[0]=='LWS': 
        #     data = sle.height_to_EWH(data).reshape(180,360)
        slf_tr = sle.run_SLE(data,name+'_tr')
        sle.plot_contour_local(X,Y,slf_tr,fname=name+'_tr',
                               save=True, savename=name+'_tr',
                               path=pwd+'plot_SLE/')
    
        data = np.array(ds.best_unc[ic_idx,:,:])
        if dataset=='IMB':
            data=-data
        # if name[ifile].split('_')[0]=='LWS':
        #     data = sle.height_to_EWH(data).reshape(180,360)
        slf_unc = sle.run_SLE(data,name+'_unc')
        sle.plot_contour_local(X,Y,np.abs(slf_unc),fname=name+'_unc', # Unc should always be positive!
                                        save=True, savename=name+'_unc',
                                        path=pwd+'plot_SLE/')
    
        
        if ifile ==0:
            df_tr = pd.DataFrame(slf_tr.flatten(),columns=[name])
            df_unc = pd.DataFrame(slf_unc.flatten(),columns=[name])
        else:
            df_tr[name]=slf_tr.flatten()
            df_unc[name]=slf_unc.flatten()
            
    #% %     
    path = '/Volumes/LaCie_NIOZ/data/barystatic/results/{}-{}/'.format(t0,t1)
    pd.to_pickle(df_tr, path+'SLF_trend_{}-{}.p'.format(t0,t1))   
    pd.to_pickle(df_unc, path+'SLF_temporal_unc_{}-{}.p'.format(t0,t1))   


