#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 16:53:55 2021

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


import matplotlib.path as mpath

#%% open datasets and merge IMBIE and JPL
path='/Volumes/LaCie_NIOZ/data/barystatic/intrinsic_unc/use/'
flist=sl.get_filelist(path,ext='*.nc')
# used IMBIE_v2 and no buf of JPL:
# flist=flist[1:3]

# ds1=xr.open_dataset(flist[1])
# print(ds1)
names=[]
for f in flist:
    ds=xr.open_dataset(f)
    # print(ds)
    n=np.array(ds.name)
    for nn in n:
        print(nn)
        names.append(nn)


#%% combine
time=np.array(ds.time)
timey=np.array(ds.year)
SL_mm=np.zeros((len(names),len(time),180,360))
SL_mm.fill('nan')
SL_EWH=np.full_like(SL_mm,np.nan)

SL_mm_y=np.zeros((len(names),len(timey),180,360))
SL_mm_y.fill('nan')
SL_EWH_y=np.full_like(SL_mm_y,np.nan)

for ifile, f in enumerate(flist):
    ds=xr.open_dataset(f)
    if ifile==0:
        j=len(ds.name)
        SL_mm[ifile:j,:,:,:]=np.array(ds.SL_mm)
        SL_EWH[ifile:j,:,:,:]=np.array(ds.SL_EWH)
        SL_mm_y[ifile:j,:,:,:]=np.array(ds.SL_mm_y)
        SL_EWH_y[ifile:j,:,:,:]=np.array(ds.SL_EWH_y)
        
    else:
        ifile=j
        j=len(ds.name)+j
        SL_mm[ifile:j,:,:,:]=np.array(ds.SL_mm)
        SL_EWH[ifile:j,:,:,:]=np.array(ds.SL_EWH)
        SL_mm_y[ifile:j,:,:,:]=np.array(ds.SL_mm_y)
        SL_EWH_y[ifile:j,:,:,:]=np.array(ds.SL_EWH_y)

#%% save 
da=xr.Dataset(data_vars={'SL_mm':(('name','time','lat','lon'),SL_mm),
                           'SL_EWH':(('name','time','lat','lon'),SL_EWH),
                           'SL_mm_y':(('name','year','lat','lon'),SL_mm_y),
                           'SL_EWH_y':(('name','year','lat','lon'),SL_EWH_y),

                          
                         },
                           coords={'lat':ds.lat,
                                   'lon':ds.lon,
                                   'time':ds.time,
                                   'tdec':ds.tdec,
                                   'year':ds.year,
                                   'name':names})
da['SL_mm'].attrs['units']='mm of sea level'
da['SL_mm_y'].attrs['units']='mm of sea level'
da['SL_EWH'].attrs['units']='mm of Equivalent Water Thickness'
da['SL_EWH_y'].attrs['units']='mm of Equivalent Water Thickness'
da['SL_mm'].attrs['long_name']='Monthly ocean mass change in mm of sea level height'
da['SL_EWH'].attrs['long_name']='Monthly ocean mass change in mm of equivalent water thickness'
da['SL_mm_y'].attrs['long_name']='Yearly averages of ocean mass change in mm of sea level height'
da['SL_EWH_y'].attrs['long_name']='Yearly averages of ocean mass change in mm of equivalent water thickness'

da.attrs['script']='intrinsic_unc-SLE.py'
da.attrs['Comment']='different sources of barystatic contributions standardized to use as input for the Sl Equation'
da.attrs['Author']='Carolina M.L. Camargo'
da.attrs['reference']='Camargo et al. 2021'
da.attrs['date_created']=str(dt.datetime.now())

path='/Volumes/LaCie_NIOZ/data/barystatic/intrinsic_unc/use/'+'comb/'
da.to_netcdf(path+'SL_unc_JPL_IMBIE.nc')
