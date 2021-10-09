#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  8 10:07:02 2021

Instead of using the uncertainity for each region, we will use the total uncertainty

@author: ccamargo
"""

import os
import xarray as xr
import numpy as np
import warnings
import sys
sys.path.append("/Users/ccamargo/Documents/py_scripts/")
import utils_SL as sl 

import utils_SLE as sle
import pandas as pd

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.gridspec import GridSpec
import cmocean as cm
from cartopy import crs as ccrs, feature as cfeature

from scipy.integrate import simps,trapz#,cumtrapz,romb

#%%# % % get ANT dataset
# IMBIE 2018 Antartic Dataset
#This spreadsheet contains the IMBIE-2018 datasets for Antarctica (worksheet 1),
#Antarctic Peninsula (worksheet 2), East Antarctica (worksheet 3) and West Antarctica (worksheet 4). 
#Each worksheet includes data on monthly cumulative ice sheet mass changes and their estimated uncertainty. 
#The data are expressed in units of mass (Gigatons – columns B and C) 
#and in units of equivalent mean global sea level rise (millimetres – columns D and E).

df=pd.read_excel('/Volumes/LaCie_NIOZ/data/barystatic/original/ANT_IMBIE.xlsx',
                  #sheet_name=[0,1] # load first and second sheet as a dict of df
                  sheet_name=0 # sheet 1
                  )

df.columns
i,j=df.shape
data=np.zeros((i,j))
c=[None]*(j+1)
for j,col in enumerate(df.columns):
    data[:,j]=np.array(df[col])
    c[j]=col

ant=data[:,4] # uncertainity
years=np.array(df.Year)
tstring=sle.from_tdec_to_tstring(years)

#% % Now get data for each region:
sheets=[1,2,3]
icesheets=['AP','EA','WA']
ant_is=np.zeros((len(sheets),len(ant)))
glb=[0,0,0]
for l,n in enumerate(sheets):
    print(n)
    df=pd.read_excel('/Volumes/LaCie_NIOZ/data/barystatic/original/ANT_IMBIE.xlsx',
                      #sheet_name=[0,1] # load first and second sheet as a dict of df
                      sheet_name=n # sheet 1
                      )
    df.columns
    i,j=df.shape
    data=np.zeros((i,j))
    c=[None]*(j+1)
    for j,col in enumerate(df.columns):
        data[:,j]=np.array(df[col])
        c[j]=col
    
    ant_is[l,:]=data[:,4] # sea level contr uncertianity (mm)


#%%  make regional time series
ice_in=np.array(ant_is)
path_to_mask='/Volumes/LaCie_NIOZ/data/barystatic/masks/barystatic_mask2.nc'
# path_to_mask='/Volumes/LaCie_NIOZ/data/barystatic/masks/rignot_basins/ANT_Basins_IMBIE2_v1.6/ANT_Basins_IMBIE2_v1.6.nc'
# path_to_mask='/Volumes/LaCie_NIOZ/PhD/Barystatic/antarctica/rignot2019/basins/AIS_basins_v2.nc'
path_to_save='/Volumes/LaCie_NIOZ/data/barystatic/intrinsic_unc/',

#% %
#% % open ANT mask:
ds=xr.open_dataset(path_to_mask)
ds
# on IMBIE we don't have island separation, so we will use 'mask_noisland':
mask=np.array(ds.mask10) # using zwally_IceSheets_1x1_etopo_cor.nc
print(ds.code_mask10)
# we want codes 10-12
mask[np.where(mask>12)]='nan'
mask[np.where(mask<10)]='nan'
np.unique(mask)
grid=np.zeros((mask.shape))
#% %
# mask=np.array(ds.mask_noisland)
mask=np.hstack(mask)
lat=np.array(ds.lat)
lon=np.array(ds.lon)
llon,llat=np.meshgrid(lon,lat)
#llat,llon=np.meshgrid(lat,lon) THIS IS WRONG!!!!!!!!!!!!!!!
llat=np.hstack(llat)
llon=np.hstack(llon)

df=pd.DataFrame(llon,columns=['lon'])
df['lat']=llat
df['mask']=mask

#%% get grid area
grid_area=sl.get_grid_area(grid)
df['area']=np.hstack(grid_area)
area_per_grid=np.hstack(grid_area)
#%% test for one time step
k=1
ice_slc=np.zeros((len(llat)))
    
# AIS
mean_load=np.array(ant[k])
region_area=df['area'].where(np.isfinite(df['mask'])).sum()
load_per_area=np.array(mean_load/region_area)
ice_slc[np.where(np.isfinite(mask))]=load_per_area*area_per_grid[np.where(np.isfinite(mask))]

print(np.sum(ant[k]))
print(np.sum(ice_slc))

#%% loop over time
# use mean value of AIS
# codes=[2,4,3]
data_reg=np.zeros((len(years),len(ds.lat),len(ds.lon)))

# Loop over time:
for k in range(0,len(years)):
    print(k)
    #Divide global number by total area of each region to obtain the mean number:
    ice_slc=np.zeros((len(llat)))
    
    # AIS
    mean_load=np.array(ant[k])
    region_area=df['area'].where(np.isfinite(df['mask'])).sum()
    load_per_area=np.array(mean_load/region_area)
    ice_slc[np.where(np.isfinite(mask))]=load_per_area*area_per_grid[np.where(np.isfinite(mask))]
    
    
    data_reg[k,:,:]=np.array(ice_slc).reshape((len(ds.lat),len(ds.lon))) 
df['ice_slc']=ice_slc   

#%% make dataset

da=xr.Dataset(data_vars={'mask':(('lat','lon'),mask.reshape(180,360)),
                          'ant_slc_reg':(('time','lat','lon'),-data_reg),
                          'ant_glb':(('time'),-ant)
                         
                        },
                          coords={'lat':ds.lat,
                                  'lon':ds.lon,

                                  'time':years,
                                  'reg':['AIS']
                                  })


#% %
da.attrs['mask_exp']='Mask for Antarctica Basisns: 0 = land + ocean; 3 = Ant Peninsula; 1 = West Ant; 2 = East Ant'
da.attrs['Metadata0']='ant_glb: global values of cumlative sea-level change for entire Antarctica (sheet 0 of IMBIE dataset)'
da.attrs['Metadata1']='ant_slc_reg: Sea-level change UNCERTAINTY from Regional Cumulative sea level contribution (mm), done by redistributing the global values of cumulative slc'
da.attrs['source'] = 'IMBIE dataset 2018 (imbie_dataset-2018_07_23).'

da.attrs['script']='intrinsic_unc-IMBIE-prep_v2'


# da.to_netcdf('/Users/ccamargo/Documents/PhD/Barystatic/imbie/imbie_ant_reg.nc')
da.to_netcdf(path_to_save[0]+'AIS_IMB.nc')

#%% GRE
#%% Greenland

# # IMBIE 2019 Greendland Dataset
# #This spreadsheet contains the IMBIE-2019 datasets for Greenland, 
# #which includes data on the annual rate of change and cumulative change in Greenland’s ice sheet mass, 
# #its surface mass balance and ice discharge anomalies, and their estimated uncertainty. 

# # Sheet 2: equivalent mean global sea level rise (mm/yr) and sea-level rise (mm)
# # and in units of equivalent mean global sea level rise (millimetres per year – sheet 2, columns B, C, F, G, J and K, 
# #and millimetres – sheet 2, columns D, E, H, I, L and M).
df=pd.read_excel('/Volumes/LaCie_NIOZ/data/barystatic/original/GRE_IMBIE.xlsx',
                  #sheet_name=[0,1] # load first and second sheet as a dict of df
                  sheet_name=1 # sheet 2
                  )
i,j=df.shape
data=np.zeros((i,j))
c=[None]*(j+1)
for j,col in enumerate(df.columns):
    data[:,j]=np.array(df[col])
    c[j]=col

years=np.array(df.Year)
print(years[144])

# From 1992 onwards UNCERTAINTIES

gre_tot=data[144:len(years),4] #Cumulative ice sheet mass change (mm sea level) 

years=years[144:len(years)]
#% % 
plt.plot(years,gre_tot,'green',label='total')

#%% Now using the Greenland Mask
#open regional grid mask:
# ds=xr.open_dataset('/Volumes/LaCie_NIOZ/PhD/Barystatic/greenland/mouginout2019/GrIS_division/GIS_basins.nc')
# ds=xr.open_dataset('/Volumes/LaCie_NIOZ/data/barystatic/masks/Mouginot_GrIS_division/Greenland_Basins_PS_v1.4.2.shp.nc')
ds=xr.open_dataset(path_to_mask)
ds
# 0:NW, 1:CE, 2:CW, 3:SW, 4:NO, 5:NE, 6:SE
ds.mask10.plot()
maskgre=np.array(ds.mask10)
# we want from 2 to 8
maskgre[np.where(maskgre>8)]='nan'
maskgre[np.where(maskgre<2)]='nan'
# remove 2 so we have it from 0-6
maskgre=maskgre-2
np.unique(maskgre[np.isfinite(maskgre)])

lat=np.array(ds.lat)
lon=np.array(ds.lon)
llon,llat=np.meshgrid(lon,lat)
#llat,llon=np.meshgrid(lat,lon) THIS IS WRONG!!!!!!!!!!!!!!!
llat=np.hstack(llat)
llon=np.hstack(llon)

df=pd.DataFrame(llon,columns=['lon'])
df['lat']=llat
df['mask']=np.hstack(maskgre)

df['area']=np.hstack(grid_area)

area_per_grid=np.hstack(grid_area)
#%% GRE TOT
data_reg = np.zeros((len(gre_tot),len(lat),len(lon)))


for k in range(0,len(years)):
    print(k)
    
    # GIS
    ice_slc=np.zeros((len(llat)))
    mean_load=np.array(gre_tot[k])
    #Divide global number by total area of each region to obtain the mean number:
    region_area=df['area'].where(np.isfinite(df['mask'])).sum()
    load_per_area=np.array(mean_load/region_area)
    ice_slc[np.isfinite(df['mask'])]=load_per_area*area_per_grid[np.isfinite(df['mask'])]
    

    #% %
    data_reg[k,:,:]=np.array(ice_slc).reshape((len(ds.lat),len(ds.lon)))
  
#%% 

#% %
da=xr.Dataset(data_vars={'gre_slc_tot':(('time','lat','lon'),-data_reg),
                         'gre_glb':(('time'),-gre_tot),
                         'mask':(('lat','lon'),maskgre),
#                         'maskdyn':(('lat','lon'),maskdyn)  ,
                         
                         },
                         coords={'lat':lat,
                                 'lon':lon,
                                 'time':years,
                                 'reg':['TOT']
                                 })



#%%
da.attrs['maskgre_exp']='Greenland mask: 0:NW, 1:CE, 2:CW, 3:SW, 4:NO, 5:NE, 6:SE Greenland Drainage division from Mouginot and Rignot, 2019'
da.attrs['Metadata']='gre_glb: global values of cumlative sea-level change UNCERTAINTY for entire ice sheet (column E of excel)'
da.attrs['units']='mm'
da.attrs['regional']='Sea-level change UNC from Regional Cumulative sea level contribution (mm), done by redistributing the global values of cumulative slc'

da.attrs['source'] = 'IMBIE dataset 2019 (imbie_dataset_greenland_dynamics-2020_02_28).'
da.attrs['script']='intrinsic_unc-IMBIE-prep_v2'
# da.to_netcdf('/Users/ccamargo/Documents/PhD/Barystatic/imbie/imbie_gre_reg.nc')
# da.to_netcdf('/Volumes/LaCie_NIOZ/PhD/Barystatic/imbie/imbie_gre_reg.nc')
da.to_netcdf(path_to_save[0]+'GIS_IMB.nc')


#%% Standardize
#%%
import datetime as dt
import utils_SLE_v2 as sle
#%%
#%% Description
# from the barystatic_test_units.py we saw that:
# datasets from IMBIE, Rignot 2019, MOUGINOT 2019 are in mm of SL

flist=sl.get_filelist(path_to_save[0],ext='*.nc')
flist=[f for f in flist if f.split('.')[-2].split('_')[-1]=='IMB']

print(flist)

names=['AIS_IMB',
       'GIS_IMB'] 
variables = ['ant_slc_reg','gre_slc_tot']
#%%
# time range: months 1993-2020.08
time,time2=sl.make_date('01-01-1993',(27*12)+8)
tdec,tdec0=sl.get_dec_time(time2,ns=1e-6)
t=[dt.datetime.utcfromtimestamp(t.astype(int) * 1e-6) for t in time2]
ty=np.array([dt.datetime.utcfromtimestamp(t.astype(int) * 1e-6).timetuple().tm_year
                                for t in time2])
tm=np.array([dt.datetime.utcfromtimestamp(t.astype(int) * 1e-6).timetuple().tm_mon
                                for t in time2])

# yearly
timey=np.arange(1993,2021)
idx=np.zeros((len(timey)))
for i,year in enumerate(timey):
    idx[i]=np.where(ty==year)[0][0]
    
SL_mm=np.zeros((len(names),len(time),180,360))
SL_mm.fill('nan')
SL_EWH=np.full_like(SL_mm,np.nan)

SL_mm_y=np.zeros((len(names),len(timey),180,360))
SL_mm_y.fill('nan')
SL_EWH_y=np.full_like(SL_mm_y,np.nan)

for iname,f in enumerate(flist):
    ds=xr.open_dataset(f)
    tdec_local=np.array(ds.time)
    time_local,time2_local=sl.make_date('01-01-'+str(int(tdec_local[0])),len(tdec_local))
    ds['time']=time_local
    ds_y=ds.groupby('time.year').mean(dim='time')
    
    name=variables[iname]
    lon=np.array(ds.lon)
    lat=np.array(ds.lat)
    timey_local=np.array(ds_y.year)
    
    # find index for yearly:
    idx_local=np.zeros((len(timey)))
    idx_local.fill('nan')
    # idx_local=np.full_like(timey,np.nan)
    for i,year in enumerate(timey):
        if year<=np.max(timey_local) and \
            year>=np.min(timey_local):
            idx_local[i]=np.where(timey_local==year)[0][0]
    
    idx_0y=np.array(np.where(np.isfinite(idx_local)))[0,0]
    idx_00y=int(idx_local[np.isfinite(idx_local)][0])
    idx_1y=np.array(np.where(np.isfinite(idx_local)))[0,-1]+1
    
    
    # find index for monthly:
    ty_local=np.array([dt.datetime.utcfromtimestamp(t.astype(int) * 1e-6).timetuple().tm_year
                                    for t in time2_local])
    tm_local=np.array([dt.datetime.utcfromtimestamp(t.astype(int) * 1e-6).timetuple().tm_mon
                                    for t in time2_local])
    
    idx_local_m=np.zeros((len(ty)))
    idx_local_m.fill('nan')
    t2_local=[dt.datetime.utcfromtimestamp(t.astype(int) * 1e-6) for t in time2_local]
    for i,iy in enumerate(t):
        year=iy.year
        if np.any(ty_local==year):
            month=iy.month
            if np.any(tm_local==month):
                # print(year)
                # print(month)
                jdx=np.array(np.where(year==ty_local))[0]
                jdx2=np.array(np.where(month==tm_local))[0]
                for jd in jdx2:
                    if np.any(jdx==jd):
                        j=jdx[jdx==jd]
                        idx_local_m[i]=j
    
    # Get data
    data=np.array(ds[name]) # mmm of SL height
    datay=np.array(ds_y[name]) 
    # Put data in main array
    SL_mm_y[iname,idx_0y:idx_1y,:,:]=datay[idx_00y:len(timey_local),:,:]
    for i,j in enumerate(idx_local_m):
        if np.isfinite(j):
            j=int(j)
            SL_mm[iname,i,:,:]=data[j,:,:]
    # Transform data: 
    for i in range(len(data)):
        data[i,:,:]=sle.height_to_EWH(data[i,:,:]).reshape(180,360)
    for i in range(len(datay)):
        datay[i,:,:]=sle.height_to_EWH(datay[i,:,:]).reshape(180,360)
    # Put transformed data in main array:
    SL_EWH_y[iname,idx_0y:idx_1y,:,:]=datay[idx_00y:len(timey_local),:,:]
    for i,j in enumerate(idx_local_m):
        if np.isfinite(j):
            j=int(j)
            SL_EWH[iname,i,:,:]=data[j,:,:]


#%%
da=xr.Dataset(data_vars={'SL_mm':(('name','time','lat','lon'),SL_mm),
                              'SL_EWH':(('name','time','lat','lon'),SL_EWH),
                              'SL_mm_y':(('name','year','lat','lon'),SL_mm_y),
                              'SL_EWH_y':(('name','year','lat','lon'),SL_EWH_y),
                              'mask_gre':(('lat','lon'),maskgre),
                              'mask_ant':(('lat','lon'),mask.reshape(180,360)),
                             
                            },
                              coords={'lat':ds.lat,
                                      'lon':ds.lon,
                                      'time':time2,
                                      'tdec':tdec,
                                      'year':timey,
                                      'name':names})
da['SL_mm'].attrs['units']='mm of sea level'
da['SL_mm_y'].attrs['units']='mm of sea level'
da['SL_EWH'].attrs['units']='mm of Equivalent Water Thickness'
da['SL_EWH_y'].attrs['units']='mm of Equivalent Water Thickness'
da['SL_mm'].attrs['long_name']='Monthly ocean mass change in mm of sea level height'
da['SL_EWH'].attrs['long_name']='Monthly ocean mass change in mm of equivalent water thickness'
da['SL_mm_y'].attrs['long_name']='Yearly averages of ocean mass change in mm of sea level height'
da['SL_EWH_y'].attrs['long_name']='Yearly averages of ocean mass change in mm of equivalent water thickness'

da['name'].attrs['long_name']='Name of the source region/contribution of ocean mass change'
#% %
da.attrs['metadata']='Ocean mass sea level changes UNCERTANTIES (in mm and EWH), from IMBIE ANT and GRE'
da.attrs['sources']='IMBIE Team, 2020; 2018, Nature '

da.attrs['script']='intrinsic_unc-IMBIE-prep.py'
da.attrs['Comment']='different sources of barystatic contributions standardized to use as input for the Sl Equation'
da.attrs['Author']='Carolina M.L. Camargo'
da.attrs['reference']='Camargo et al. 2021'
da.attrs['date_created']=str(dt.datetime.now())

path='/Volumes/LaCie_NIOZ/data/barystatic/intrinsic_unc/use/'
da.to_netcdf(path+'IMBIE_v2.nc')

