#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Thu Mar 11 17:29:48 2021

Combine different datasets

Make units and time period standard

@author: ccamargo
"""


import numpy as np
import xarray as xr
import sys
sys.path.append("/Users/ccamargo/Documents/py_scripts/")
import utils_SL as sl 
import utils_SLE_v2 as sle

# from netCDF4 import Dataset

import datetime as dt

#%% Description
# from the barystatic_test_units.py we saw that:
# datasets from IMBIE, Rignot 2019, MOUGINOT 2019 are in mm of SL
# dataset from Zemp 2019 is in mm of SL 
# dataset from PCR-GLOBWB is in dcm of Equivalent Water Thickness [x100]
# dataset from WaterGAP is in mm of Equivalent Water Thickness (kg m-2)
# datasets from CSR and JPL Mascons are in cm of Equivalent Water Thickness [x10]
path='/Volumes/LaCie_NIOZ/data/barystatic/use/'

flist=sl.get_filelist(path,ext='*.nc')
print(flist)

names=['AIS_IMB','AIS_R19','AIS_R19_basins',
       'GLWS_ZMP',
       'AIS_300_CSR','GIS_300_CSR',
       'LWS_CSR','GLWS_CSR',
       'TCWS_CSR',
       'AIS_CSR','GIS_CSR',
       'AIS_300_JPL','GIS_300_JPL',
       'LWS_JPL','GLWS_JPL','TCWS_JPL',
       'AIS_JPL','GIS_JPL',
       'AIS_proj_CSR','GIS_proj_CSR',
       'AIS_proj_JPL','GIS_proj_JPL',
       'GIS_IMB','GIS_M19',
       'LWS_GWB',
       'LWS_WGP_gl','GLWS_WGP_gl',
       'TCWS_WaterGAP'] 

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

#%% 1. ANT_IMB
ifile=0
iname=0
f=flist[ifile]
name=names[iname]
print(name)
ds=xr.open_dataset(f)
name='ant_slc_reg'
tdec_local=np.array(ds.time)
time_local,time2_local=sl.make_date('01-01-'+str(int(tdec_local[0])),len(tdec_local))
da=xr.Dataset(data_vars={name:(('time','lat','lon'),ds[name]),
                        },
                         coords={'lat':ds.lat,
                                 'lon':ds.lon,
                                 'time':time_local})
ds_y=da.groupby('time.year').mean(dim='time')

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
data=np.array(da[name]) # mmm of SL height
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
        
#%% 2. ANT_R19
ifile=ifile+1
iname=iname+1
f=flist[ifile]
name=names[iname]
print(name)
ds=xr.open_dataset(f) # yearly
name='MB'

time=np.array(ds.time)
lon=np.array(ds.lon)
lat=np.array(ds.lat)
data=np.array(ds[name])
datay=np.array(ds[name])
# Find local index
idx_local=np.zeros((len(timey)))
idx_local.fill('nan')
for i,year in enumerate(timey):
    if year<=np.max(time):
        idx_local[i]=np.where(time==year)[0][0]

idx_0=int(idx_local[0])
idx_1=np.array(np.where(np.isfinite(idx_local)))[0,-1]+1

# Put data in main array
SL_mm_y[iname,0:idx_1,:,:]=data[idx_0:len(time),:,:]
for i,j in zip(idx,idx_local):
    if np.isfinite(j):
        i=int(i);j=int(j)
        SL_mm[iname,i,:,:]=data[j,:,:]

# Transform data: 
for i in range(len(data)):
    data[i,:,:]=sle.height_to_EWH(data[i,:,:]).reshape(180,360)
for i in range(len(datay)):
    datay[i,:,:]=sle.height_to_EWH(datay[i,:,:]).reshape(180,360)

# Put transformed data in main array
SL_EWH_y[iname,0:idx_1,:,:]=data[idx_0:len(time),:,:]

for i,j in zip(idx,idx_local):
    if np.isfinite(j):
        i=int(i);j=int(j)
        SL_EWH[iname,i,:,:]=data[j,:,:]
#%% 2. ANT_R19
ifile=ifile+1
iname=iname+1
f=flist[ifile]
name=names[iname]
print(name)
ds=xr.open_dataset(f) # yearly
name='MB'

if np.nanmin(ds.lon)<0:
    lon=sl.from_180_to_360(np.array(ds.lon))
    ds = ds.assign_coords(lon=lon)
    ds = ds.sortby('lon')

time=np.array(ds.time)
lon=np.array(ds.lon)
lat=np.array(ds.lat)

data=np.array(ds[name])
datay=np.array(ds[name])
# Find local index
idx_local=np.zeros((len(timey)))
idx_local.fill('nan')
for i,year in enumerate(timey):
    if year<=np.max(time):
        idx_local[i]=np.where(time==year)[0][0]

idx_0=int(idx_local[0])
idx_1=np.array(np.where(np.isfinite(idx_local)))[0,-1]+1

# Put data in main array
SL_mm_y[iname,0:idx_1,:,:]=data[idx_0:len(time),:,:]
for i,j in zip(idx,idx_local):
    if np.isfinite(j):
        i=int(i);j=int(j)
        SL_mm[iname,i,:,:]=data[j,:,:]

# Transform data: 
for i in range(len(data)):
    data[i,:,:]=sle.height_to_EWH(data[i,:,:]).reshape(180,360)
for i in range(len(datay)):
    datay[i,:,:]=sle.height_to_EWH(datay[i,:,:]).reshape(180,360)

# Put transformed data in main array
SL_EWH_y[iname,0:idx_1,:,:]=data[idx_0:len(time),:,:]

for i,j in zip(idx,idx_local):
    if np.isfinite(j):
        i=int(i);j=int(j)
        SL_EWH[iname,i,:,:]=data[j,:,:]

#%% 3. GLA ZMP
ifile=ifile+1
iname=iname+1
f=flist[ifile]
name=names[iname]
print(name)
ds=xr.open_dataset(f) # already yearly
# ds_y=ds.groupby('time.year').mean(dim='time')

name='GMC_AW_01' # no adjacent glaciers to GRE and ANT
time=np.array(ds.time)
lon=np.array(ds.lon)
lat=np.array(ds.lat)
data=np.array(ds[name]) #%% mm of SL height
datay=np.array(ds[name]) # data is already in years

# Find local index
idx_local=np.zeros((len(timey)))
idx_local.fill('nan')
for i,year in enumerate(timey):
    if year<=np.max(time):
        idx_local[i]=np.where(time==year)[0][0]

idx_0=int(idx_local[0])
idx_1=np.array(np.where(np.isfinite(idx_local)))[0,-1]+1

# Put data in main array
SL_mm_y[iname,0:idx_1,:,:]=data[idx_0:len(time),:,:]
for i,j in zip(idx,idx_local):
    if np.isfinite(j):
        i=int(i);j=int(j)
        SL_mm[iname,i,:,:]=data[j,:,:]

# Transform data: 
for i in range(len(data)):
    data[i,:,:]=sle.height_to_EWH(data[i,:,:]).reshape(180,360)
for i in range(len(datay)):
    datay[i,:,:]=sle.height_to_EWH(datay[i,:,:]).reshape(180,360)

# Put transformed data in main array
SL_EWH_y[iname,0:idx_1,:,:]=data[idx_0:len(time),:,:]

for i,j in zip(idx,idx_local):
    if np.isfinite(j):
        i=int(i);j=int(j)
        SL_EWH[iname,i,:,:]=data[j,:,:]
#%% 4. CSR-300km
ifile=ifile+1
f=flist[ifile]
ds=xr.open_dataset(f)
ds_y=ds.groupby('time.year').mean(dim='time')

# % % find time index
time_local=np.array(ds.time)
tdec_local,tdec0=sl.get_dec_time(time_local,ns=1e-9)
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
ty_local=np.array([dt.datetime.utcfromtimestamp(t.astype(int) * 1e-9).timetuple().tm_year
                                for t in time_local])
tm_local=np.array([dt.datetime.utcfromtimestamp(t.astype(int) * 1e-9).timetuple().tm_mon
                                for t in time_local])

idx_local_m=np.zeros((len(ty)))
idx_local_m.fill('nan')
t2_local=[dt.datetime.utcfromtimestamp(t.astype(int) * 1e-9 ) for t in time_local]
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
i=225
idx_local_m[i+1]=idx_local_m[i]
idx_local_m[i]=idx_local_m[i]-1
i=267
idx_local_m[i+1]=idx_local_m[i]
idx_local_m[i]=idx_local_m[i]-1

## ANT_CSR 300km
iname=iname+1
name=names[iname]
print(name)
name='ant'
 # Water thickness
data=np.array(ds[name])*10 # mm of water thickness
datay=np.array(ds_y[name])*10 
# data[np.where(data==-9990)]=np.nan
# datay[np.where(datay==-9990)]=np.nan
# Put data in main array
SL_EWH_y[iname,idx_0y:idx_1y,:,:]=datay[idx_00y:len(timey_local),:,:]
for i,j in enumerate(idx_local_m):
    if np.isfinite(j):
        j=int(j)
        SL_EWH[iname,i,:,:]=data[j,:,:]
# Transform data: 
for i in range(len(data)):
    data[i,:,:]=sle.EWH_to_height(data[i,:,:]).reshape(180,360)
for i in range(len(datay)):
    datay[i,:,:]=sle.EWH_to_height(datay[i,:,:]).reshape(180,360)
# Put transformed data in main array:
SL_mm_y[iname,idx_0y:idx_1y,:,:]=datay[idx_00y:len(timey_local),:,:]
for i,j in enumerate(idx_local_m):
    if np.isfinite(j):
        j=int(j)
        SL_mm[iname,i,:,:]=data[j,:,:]
        
## GRE_CSR 300km
iname=iname+1
name=names[iname]
print(name)

name='gre'
 # Water thickness
data=np.array(ds[name])*10 # mm of water thickness
datay=np.array(ds_y[name])*10 
# data[np.where(data==-9990)]=np.nan
# datay[np.where(datay==-9990)]=np.nan
# Put data in main array
SL_EWH_y[iname,idx_0y:idx_1y,:,:]=datay[idx_00y:len(timey_local),:,:]
for i,j in enumerate(idx_local_m):
    if np.isfinite(j):
        j=int(j)
        SL_EWH[iname,i,:,:]=data[j,:,:]
# Transform data: 
for i in range(len(data)):
    data[i,:,:]=sle.EWH_to_height(data[i,:,:]).reshape(180,360)
for i in range(len(datay)):
    datay[i,:,:]=sle.EWH_to_height(datay[i,:,:]).reshape(180,360)
# Put transformed data in main array:
SL_mm_y[iname,idx_0y:idx_1y,:,:]=datay[idx_00y:len(timey_local),:,:]
for i,j in enumerate(idx_local_m):
    if np.isfinite(j):
        j=int(j)
        SL_mm[iname,i,:,:]=data[j,:,:]

## TWS_CSR
iname=iname+1
name=names[iname]
print(name)

name='tws'
 # Water thickness
data=np.array(ds[name])*10 # mm of water thickness
datay=np.array(ds_y[name])*10 
# data[np.where(data==-9990)]=np.nan
# datay[np.where(datay==-9990)]=np.nan
# Put data in main array
SL_EWH_y[iname,idx_0y:idx_1y,:,:]=datay[idx_00y:len(timey_local),:,:]
for i,j in enumerate(idx_local_m):
    if np.isfinite(j):
        j=int(j)
        SL_EWH[iname,i,:,:]=data[j,:,:]
# Transform data: 
for i in range(len(data)):
    data[i,:,:]=sle.EWH_to_height(data[i,:,:]).reshape(180,360)
for i in range(len(datay)):
    datay[i,:,:]=sle.EWH_to_height(datay[i,:,:]).reshape(180,360)
# Put transformed data in main array:
SL_mm_y[iname,idx_0y:idx_1y,:,:]=datay[idx_00y:len(timey_local),:,:]
for i,j in enumerate(idx_local_m):
    if np.isfinite(j):
        j=int(j)
        SL_mm[iname,i,:,:]=data[j,:,:]
## GLA CSR
iname=iname+1
name=names[iname]
print(name)

name='gla'
 # Water thickness
data=np.array(ds[name])*10 # mm of water thickness
datay=np.array(ds_y[name])*10 
# data[np.where(data==-9990)]=np.nan
# datay[np.where(datay==-9990)]=np.nan
# Put data in main array
SL_EWH_y[iname,idx_0y:idx_1y,:,:]=datay[idx_00y:len(timey_local),:,:]
for i,j in enumerate(idx_local_m):
    if np.isfinite(j):
        j=int(j)
        SL_EWH[iname,i,:,:]=data[j,:,:]
# Transform data: 
for i in range(len(data)):
    data[i,:,:]=sle.EWH_to_height(data[i,:,:]).reshape(180,360)
for i in range(len(datay)):
    datay[i,:,:]=sle.EWH_to_height(datay[i,:,:]).reshape(180,360)
# Put transformed data in main array:
SL_mm_y[iname,idx_0y:idx_1y,:,:]=datay[idx_00y:len(timey_local),:,:]
for i,j in enumerate(idx_local_m):
    if np.isfinite(j):
        j=int(j)
        SL_mm[iname,i,:,:]=data[j,:,:]
## BAR CSR
iname=iname+1
name=names[iname]
print(name)

name='bary_total'
 # Water thickness
data=np.array(ds[name])*10 # mm of water thickness
datay=np.array(ds_y[name])*10 
# Put data in main array
SL_EWH_y[iname,idx_0y:idx_1y,:,:]=datay[idx_00y:len(timey_local),:,:]
for i,j in enumerate(idx_local_m):
    if np.isfinite(j):
        j=int(j)
        SL_EWH[iname,i,:,:]=data[j,:,:]
# Transform data: 
for i in range(len(data)):
    data[i,:,:]=sle.EWH_to_height(data[i,:,:]).reshape(180,360)
for i in range(len(datay)):
    datay[i,:,:]=sle.EWH_to_height(datay[i,:,:]).reshape(180,360)
# Put transformed data in main array:
SL_mm_y[iname,idx_0y:idx_1y,:,:]=datay[idx_00y:len(timey_local),:,:]
for i,j in enumerate(idx_local_m):
    if np.isfinite(j):
        j=int(j)
        SL_mm[iname,i,:,:]=data[j,:,:]
#%% 5. CSR no buf
ifile=ifile+1
f=flist[ifile]
ds=xr.open_dataset(f)
ds_y=ds.groupby('time.year').mean(dim='time')
 # no need to find time index
 # same as previous
## ANT_CSR no buf
iname=iname+1
name=names[iname]
print(name)

name='ant'
 # Water thickness
data=np.array(ds[name])*10 # mm of water thickness
datay=np.array(ds_y[name])*10 
# Put data in main array
SL_EWH_y[iname,idx_0y:idx_1y,:,:]=datay[idx_00y:len(timey_local),:,:]
for i,j in enumerate(idx_local_m):
    if np.isfinite(j):
        j=int(j)
        SL_EWH[iname,i,:,:]=data[j,:,:]
# Transform data: 
for i in range(len(data)):
    data[i,:,:]=sle.EWH_to_height(data[i,:,:]).reshape(180,360)
for i in range(len(datay)):
    datay[i,:,:]=sle.EWH_to_height(datay[i,:,:]).reshape(180,360)
# Put transformed data in main array:
SL_mm_y[iname,idx_0y:idx_1y,:,:]=datay[idx_00y:len(timey_local),:,:]
for i,j in enumerate(idx_local_m):
    if np.isfinite(j):
        j=int(j)
        SL_mm[iname,i,:,:]=data[j,:,:]
        
## GRE_CSR no buf
iname=iname+1
name=names[iname]
print(name)

name='gre'
 # Water thickness
data=np.array(ds[name])*10 # mm of water thickness
datay=np.array(ds_y[name])*10 
# Put data in main array
SL_EWH_y[iname,idx_0y:idx_1y,:,:]=datay[idx_00y:len(timey_local),:,:]
for i,j in enumerate(idx_local_m):
    if np.isfinite(j):
        j=int(j)
        SL_EWH[iname,i,:,:]=data[j,:,:]
# Transform data: 
for i in range(len(data)):
    data[i,:,:]=sle.EWH_to_height(data[i,:,:]).reshape(180,360)
for i in range(len(datay)):
    datay[i,:,:]=sle.EWH_to_height(datay[i,:,:]).reshape(180,360)
# Put transformed data in main array:
SL_mm_y[iname,idx_0y:idx_1y,:,:]=datay[idx_00y:len(timey_local),:,:]
for i,j in enumerate(idx_local_m):
    if np.isfinite(j):
        j=int(j)
        SL_mm[iname,i,:,:]=data[j,:,:]

#%% 6. JPL-300km
ifile=ifile+1
f=flist[ifile]

ds=xr.open_dataset(f)
ds_y=ds.groupby('time.year').mean(dim='time')
 # no need to find time index
 # same as previous

## ANT_JPL 300km
iname=iname+1
name=names[iname]
print(name)

name='ant'
 # Water thickness
data=np.array(ds[name])*10 # mm of water thickness
datay=np.array(ds_y[name])*10 
# Put data in main array
SL_EWH_y[iname,idx_0y:idx_1y,:,:]=datay[idx_00y:len(timey_local),:,:]
for i,j in enumerate(idx_local_m):
    if np.isfinite(j):
        j=int(j)
        SL_EWH[iname,i,:,:]=data[j,:,:]
# Transform data: 
for i in range(len(data)):
    data[i,:,:]=sle.EWH_to_height(data[i,:,:]).reshape(180,360)
for i in range(len(datay)):
    datay[i,:,:]=sle.EWH_to_height(datay[i,:,:]).reshape(180,360)
# Put transformed data in main array:
SL_mm_y[iname,idx_0y:idx_1y,:,:]=datay[idx_00y:len(timey_local),:,:]
for i,j in enumerate(idx_local_m):
    if np.isfinite(j):
        j=int(j)
        SL_mm[iname,i,:,:]=data[j,:,:]
        
## GRE_JPL 300km
iname=iname+1
name=names[iname]
print(name)

name='gre'
 # Water thickness
data=np.array(ds[name])*10 # mm of water thickness
datay=np.array(ds_y[name])*10 
# Put data in main array
SL_EWH_y[iname,idx_0y:idx_1y,:,:]=datay[idx_00y:len(timey_local),:,:]
for i,j in enumerate(idx_local_m):
    if np.isfinite(j):
        j=int(j)
        SL_EWH[iname,i,:,:]=data[j,:,:]
# Transform data: 
for i in range(len(data)):
    data[i,:,:]=sle.EWH_to_height(data[i,:,:]).reshape(180,360)
for i in range(len(datay)):
    datay[i,:,:]=sle.EWH_to_height(datay[i,:,:]).reshape(180,360)
# Put transformed data in main array:
SL_mm_y[iname,idx_0y:idx_1y,:,:]=datay[idx_00y:len(timey_local),:,:]
for i,j in enumerate(idx_local_m):
    if np.isfinite(j):
        j=int(j)
        SL_mm[iname,i,:,:]=data[j,:,:]
## TWS_JPL
iname=iname+1
name=names[iname]
print(name)

name='tws'
 # Water thickness
data=np.array(ds[name])*10 # mm of water thickness
datay=np.array(ds_y[name])*10 
# Put data in main array
SL_EWH_y[iname,idx_0y:idx_1y,:,:]=datay[idx_00y:len(timey_local),:,:]
for i,j in enumerate(idx_local_m):
    if np.isfinite(j):
        j=int(j)
        SL_EWH[iname,i,:,:]=data[j,:,:]
# Transform data: 
for i in range(len(data)):
    data[i,:,:]=sle.EWH_to_height(data[i,:,:]).reshape(180,360)
for i in range(len(datay)):
    datay[i,:,:]=sle.EWH_to_height(datay[i,:,:]).reshape(180,360)
# Put transformed data in main array:
SL_mm_y[iname,idx_0y:idx_1y,:,:]=datay[idx_00y:len(timey_local),:,:]
for i,j in enumerate(idx_local_m):
    if np.isfinite(j):
        j=int(j)
        SL_mm[iname,i,:,:]=data[j,:,:]
## GLA JPL
iname=iname+1
name=names[iname]
print(name)

name='gla'
 # Water thickness
data=np.array(ds[name])*10 # mm of water thickness
datay=np.array(ds_y[name])*10 
# Put data in main array
SL_EWH_y[iname,idx_0y:idx_1y,:,:]=datay[idx_00y:len(timey_local),:,:]
for i,j in enumerate(idx_local_m):
    if np.isfinite(j):
        j=int(j)
        SL_EWH[iname,i,:,:]=data[j,:,:]
# Transform data: 
for i in range(len(data)):
    data[i,:,:]=sle.EWH_to_height(data[i,:,:]).reshape(180,360)
for i in range(len(datay)):
    datay[i,:,:]=sle.EWH_to_height(datay[i,:,:]).reshape(180,360)
# Put transformed data in main array:
SL_mm_y[iname,idx_0y:idx_1y,:,:]=datay[idx_00y:len(timey_local),:,:]
for i,j in enumerate(idx_local_m):
    if np.isfinite(j):
        j=int(j)
        SL_mm[iname,i,:,:]=data[j,:,:]
        
## BAR JPL
iname=iname+1
name=names[iname]
name='bary_total'
 # Water thickness
data=np.array(ds[name])*10 # mm of water thickness
datay=np.array(ds_y[name])*10 
# Put data in main array
SL_EWH_y[iname,idx_0y:idx_1y,:,:]=datay[idx_00y:len(timey_local),:,:]
for i,j in enumerate(idx_local_m):
    if np.isfinite(j):
        j=int(j)
        SL_EWH[iname,i,:,:]=data[j,:,:]
# Transform data: 
for i in range(len(data)):
    data[i,:,:]=sle.EWH_to_height(data[i,:,:]).reshape(180,360)
for i in range(len(datay)):
    datay[i,:,:]=sle.EWH_to_height(datay[i,:,:]).reshape(180,360)
# Put transformed data in main array:
SL_mm_y[iname,idx_0y:idx_1y,:,:]=datay[idx_00y:len(timey_local),:,:]
for i,j in enumerate(idx_local_m):
    if np.isfinite(j):
        j=int(j)
        SL_mm[iname,i,:,:]=data[j,:,:]
        
        #%% 7. JPL no buf
ifile=ifile+1
f=flist[ifile]
ds=xr.open_dataset(f)
ds_y=ds.groupby('time.year').mean(dim='time')
 # no need to find time index
 # same as previous

## ANT_JPL no buf
iname=iname+1
name=names[iname]
print(name)

name='ant'
 # Water thickness
data=np.array(ds[name])*10 # mm of water thickness
datay=np.array(ds_y[name])*10 
# Put data in main array
SL_EWH_y[iname,idx_0y:idx_1y,:,:]=datay[idx_00y:len(timey_local),:,:]
for i,j in enumerate(idx_local_m):
    if np.isfinite(j):
        j=int(j)
        SL_EWH[iname,i,:,:]=data[j,:,:]
# Transform data: 
for i in range(len(data)):
    data[i,:,:]=sle.EWH_to_height(data[i,:,:]).reshape(180,360)
for i in range(len(datay)):
    datay[i,:,:]=sle.EWH_to_height(datay[i,:,:]).reshape(180,360)
# Put transformed data in main array:
SL_mm_y[iname,idx_0y:idx_1y,:,:]=datay[idx_00y:len(timey_local),:,:]
for i,j in enumerate(idx_local_m):
    if np.isfinite(j):
        j=int(j)
        SL_mm[iname,i,:,:]=data[j,:,:]
        
## GRE_JPL no buf
iname=iname+1
name=names[iname]
print(name)

name='gre'
 # Water thickness
data=np.array(ds[name])*10 # mm of water thickness
datay=np.array(ds_y[name])*10 
# Put data in main array
SL_EWH_y[iname,idx_0y:idx_1y,:,:]=datay[idx_00y:len(timey_local),:,:]
for i,j in enumerate(idx_local_m):
    if np.isfinite(j):
        j=int(j)
        SL_EWH[iname,i,:,:]=data[j,:,:]
# Transform data: 
for i in range(len(data)):
    data[i,:,:]=sle.EWH_to_height(data[i,:,:]).reshape(180,360)
for i in range(len(datay)):
    datay[i,:,:]=sle.EWH_to_height(datay[i,:,:]).reshape(180,360)
# Put transformed data in main array:
SL_mm_y[iname,idx_0y:idx_1y,:,:]=datay[idx_00y:len(timey_local),:,:]
for i,j in enumerate(idx_local_m):
    if np.isfinite(j):
        j=int(j)
        SL_mm[iname,i,:,:]=data[j,:,:]

#%% 8. Mascons reproject
ifile=ifile+1
f=flist[ifile]
ds=xr.open_dataset(f)
ds_y=ds.groupby('time.year').mean(dim='time')
 # no need to find time index
 # same as previous
## ANT_CSR no buf
iname=iname+1
name=names[iname]
print(name)

name='AIS_CSR'
 # Water thickness
data=np.array(ds[name])*10 # mm of SL
datay=np.array(ds_y[name])*10 
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
## GRE_CSR no buf
iname=iname+1
name=names[iname]
print(name)

name='GIS_CSR'
 # Water thickness
data=np.array(ds[name])*10 # mm of SL
datay=np.array(ds_y[name])*10 
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
        
## ANT_JPL no buf
iname=iname+1
name=names[iname]
print(name)

name='AIS_JPL'
 # Water thickness
data=np.array(ds[name])*10 # mm of SL
datay=np.array(ds_y[name])*10 
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
## GRE_JPL no buf
iname=iname+1
name=names[iname]
print(name)

name='GIS_JPL'
 # Water thickness
data=np.array(ds[name])*10 # mm of SL
datay=np.array(ds_y[name])*10 
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
#%% 9. GRE imbie
ifile=ifile+1
f=flist[ifile]
ds=xr.open_dataset(f)

iname=iname+1
name=names[iname]
print(name)

name='gre_slc_tot'
tdec_local=np.array(ds.time)
time_local,time2_local=sl.make_date('01-01-'+str(int(tdec_local[0])),len(tdec_local))
da=xr.Dataset(data_vars={name:(('time','lat','lon'),ds[name]),
                        },
                         coords={'lat':ds.lat,
                                 'lon':ds.lon,
                                 'time':time_local})
ds_y=da.groupby('time.year').mean(dim='time')

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
data=np.array(da[name]) # mmm of SL height
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
#%% 10. GRE M19
ifile=ifile+1
f=flist[ifile]
ds=xr.open_dataset(f)
# ds_y=ds.groupby('time.year').mean(dim='time')
iname=iname+1
name=names[iname]
print(name)

name='SLC'
time=np.array(ds.time)
lon=np.array(ds.lon)
lat=np.array(ds.lat)
data=np.array(ds[name])
datay=np.array(ds[name]) # data is already in years
# Find local index
idx_local=np.zeros((len(timey)))
idx_local.fill('nan')
for i,year in enumerate(timey):
    if year<=np.max(time):
        idx_local[i]=np.where(time==year)[0][0]

idx_0=int(idx_local[0])
idx_1=np.array(np.where(np.isfinite(idx_local)))[0,-1]+1

# Put data in main array
SL_mm_y[iname,0:idx_1,:,:]=data[idx_0:len(time),:,:]
for i,j in zip(idx,idx_local):
    if np.isfinite(j):
        i=int(i);j=int(j)
        SL_mm[iname,i,:,:]=data[j,:,:]

# Transform data: 
for i in range(len(data)):
    data[i,:,:]=sle.height_to_EWH(data[i,:,:]).reshape(180,360)
for i in range(len(datay)):
    datay[i,:,:]=sle.height_to_EWH(datay[i,:,:]).reshape(180,360)

# Put transformed data in main array
SL_EWH_y[iname,0:idx_1,:,:]=data[idx_0:len(time),:,:]

for i,j in zip(idx,idx_local):
    if np.isfinite(j):
        i=int(i);j=int(j)
        SL_EWH[iname,i,:,:]=data[j,:,:]

#%% TWS GWB
ifile=ifile+1
f=flist[ifile]
ds=xr.open_dataset(f)
ds_y=ds.groupby('time.year').mean(dim='time')
iname=iname+1
name=names[iname]
print(name)
name='tws_no_gre_gla'

#% % find time index
time_local=np.array(ds.time)
tdec_local,tdec0=sl.get_dec_time(time_local,ns=1e-9)
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
ty_local=np.array([dt.datetime.utcfromtimestamp(t.astype(int) * 1e-9).timetuple().tm_year
                                for t in time_local])
tm_local=np.array([dt.datetime.utcfromtimestamp(t.astype(int) * 1e-9).timetuple().tm_mon
                                for t in time_local])

idx_local_m=np.zeros((len(ty)))
idx_local_m.fill('nan')
t2_local=[dt.datetime.utcfromtimestamp(t.astype(int) * 1e-9 ) for t in time_local]
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

#% % Water thickness original (dcm of water thickness)
data=np.array(ds[name])*100 # dcm->mm of water thickness
datay=np.array(ds_y[name])*100 
# Put data in main array
SL_EWH_y[iname,idx_0y:idx_1y,:,:]=datay[idx_00y:len(timey_local),:,:]
for i,j in enumerate(idx_local_m):
    if np.isfinite(j):
        j=int(j)
        SL_EWH[iname,i,:,:]=data[j,:,:]
# Transform data: 
for i in range(len(data)):
    data[i,:,:]=sle.EWH_to_height(data[i,:,:]).reshape(180,360)
for i in range(len(datay)):
    datay[i,:,:]=sle.EWH_to_height(datay[i,:,:]).reshape(180,360)
# Put transformed data in main array:
SL_mm_y[iname,idx_0y:idx_1y,:,:]=datay[idx_00y:len(timey_local),:,:]
for i,j in enumerate(idx_local_m):
    if np.isfinite(j):
        j=int(j)
        SL_mm[iname,i,:,:]=data[j,:,:]
#%% TWS WGP
ifile=ifile+1
f=flist[ifile]
ds=xr.open_dataset(f,decode_times=False)
name='tws_no_gregla'

#% % find time index
tdec_local=np.array(ds.time)
time_local,time2_local=sl.make_date('01-01-1948',len(tdec_local))
da=xr.Dataset(data_vars={name:(('time','lat','lon'),ds[name]),
                        },
                         coords={'lat':lat,
                                 'lon':lon,
                                 'time':time_local})
ds_y=da.groupby('time.year').mean(dim='time')

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

## TWS WGP
iname=iname+1
name=names[iname]
print(name)

name='tws_no_gregla'
da=xr.Dataset(data_vars={name:(('time','lat','lon'),ds[name]),
                        },
                          coords={'lat':lat,
                                  'lon':lon,
                                  'time':time_local})
ds_y=da.groupby('time.year').mean(dim='time')
data=np.array(ds[name])# mm of water thickness
datay=np.array(ds_y[name])
# Put data in main array
SL_EWH_y[iname,idx_0y:idx_1y,:,:]=datay[idx_00y:len(timey_local),:,:]
for i,j in enumerate(idx_local_m):
    if np.isfinite(j):
        j=int(j)
        SL_EWH[iname,i,:,:]=data[j,:,:]
# Transform data: 
for i in range(len(data)):
    data[i,:,:]=sle.EWH_to_height(data[i,:,:]).reshape(180,360)
for i in range(len(datay)):
    datay[i,:,:]=sle.EWH_to_height(datay[i,:,:]).reshape(180,360)
# Put transformed data in main array:
SL_mm_y[iname,idx_0y:idx_1y,:,:]=datay[idx_00y:len(timey_local),:,:]
for i,j in enumerate(idx_local_m):
    if np.isfinite(j):
        j=int(j)
        SL_mm[iname,i,:,:]=data[j,:,:]
        
## GLA WGP
iname=iname+1
name=names[iname]
print(name)

name='tws_only_gla'
da=xr.Dataset(data_vars={name:(('time','lat','lon'),ds[name]),
                        },
                         coords={'lat':lat,
                                 'lon':lon,
                                 'time':time_local})
ds_y=da.groupby('time.year').mean(dim='time')
#% % get data
data=np.array(ds[name])# mm of water thickness
datay=np.array(ds_y[name])
# Put data in main array
SL_EWH_y[iname,idx_0y:idx_1y,:,:]=datay[idx_00y:len(timey_local),:,:]
for i,j in enumerate(idx_local_m):
    if np.isfinite(j):
        j=int(j)
        SL_EWH[iname,i,:,:]=data[j,:,:]
# Transform data: 
for i in range(len(data)):
    data[i,:,:]=sle.EWH_to_height(data[i,:,:]).reshape(180,360)
for i in range(len(datay)):
    datay[i,:,:]=sle.EWH_to_height(datay[i,:,:]).reshape(180,360)
# Put transformed data in main array:
SL_mm_y[iname,idx_0y:idx_1y,:,:]=datay[idx_00y:len(timey_local),:,:]
for i,j in enumerate(idx_local_m):
    if np.isfinite(j):
        j=int(j)
        SL_mm[iname,i,:,:]=data[j,:,:]
        
##LWS WGP
iname=iname+1
name=names[iname]
print(name)

name='tws_no_gre'
da=xr.Dataset(data_vars={name:(('time','lat','lon'),ds[name]),
                        },
                         coords={'lat':lat,
                                 'lon':lon,
                                 'time':time_local})
ds_y=da.groupby('time.year').mean(dim='time')
#% % get data
data=np.array(ds[name])# mm of water thickness
datay=np.array(ds_y[name])
# Put data in main array
SL_EWH_y[iname,idx_0y:idx_1y,:,:]=datay[idx_00y:len(timey_local),:,:]
for i,j in enumerate(idx_local_m):
    if np.isfinite(j):
        j=int(j)
        SL_EWH[iname,i,:,:]=data[j,:,:]
# Transform data: 
for i in range(len(data)):
    data[i,:,:]=sle.EWH_to_height(data[i,:,:]).reshape(180,360)
for i in range(len(datay)):
    datay[i,:,:]=sle.EWH_to_height(datay[i,:,:]).reshape(180,360)
# Put transformed data in main array:
SL_mm_y[iname,idx_0y:idx_1y,:,:]=datay[idx_00y:len(timey_local),:,:]
for i,j in enumerate(idx_local_m):
    if np.isfinite(j):
        j=int(j)
        SL_mm[iname,i,:,:]=data[j,:,:]
        
#%%
#%%
#%%

#%% make and save data array
ds=xr.open_dataset('/Volumes/LaCie_NIOZ/data/barystatic/masks/barystatic_mask2.nc')
ds

#%%
da=xr.Dataset(data_vars={'SL_mm':(('name','time','lat','lon'),SL_mm),
                         'SL_EWH':(('name','time','lat','lon'),SL_EWH),
                         'SL_mm_y':(('name','year','lat','lon'),SL_mm_y),
                         'SL_EWH_y':(('name','year','lat','lon'),SL_EWH_y),
                         'mask':(('lat','lon'),ds['mask6']),
                         'mask2':(('lat','lon'),ds['mask12']),
                         
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
# da['time'].attrs['long_name']='Time'
# da['time'].attrs['standard_time']='time'
# da['time'].attrs['units']='months since january-1993'
# da['time'].attrs['calendar']='gregorian'
# da['tdec'].attrs['long_name']='Decimal time'
# da['tdec'].attrs['units']='months since january-1993'
# da['tdec'].attrs['calendar']='gregorian'
# da['year'].attrs['long_name']='Years'
# da['year'].attrs['datasets']='R19,M19 and ZMP were originally yearly data. Other datasets have been averaged to obtain an year mean.'
# da['tdec'].attrs['units']='years 1993'
da['name'].attrs['long_name']='Name of the source region/contribution of ocean mass change'
#% %
da.attrs['metadata']='Ocean mass sea level changes (in mm and EWH), from differente contributions (ANT, GRE, TWS, GLA) \
    and sources (IMBIE, Rignot 2019, Mouginot 2019, WaterGAP, PCR-GLOBWB, Zemp2019, and CSr and JPL Mascons)'
da.attrs['sources']='IMB: IMBIE Report (IMBIE Team, 2018 and 2020, Nature);\
                    R19: Rignot et al. 2019, PNAS;\
                    M19: Mouginot et al. 2019, PNAS;\
                    WGP: WaterGAP 22d with 70% irrigation Caceres et al. 2020, Hydrol. Earth. Sci\
                    GWB: PCR-GLobWB 2 from Sutanudjaja et al, 2018, Geosci. Model Dev.\
                    ZEMP: Zemp et al., 2019, Nature,\
                    CSR: CSR GRACE and GRACE-FO RL06 v02 from Save et al., 2020\
                    JPL: JPL GRACE and GRACE-FO RL06 v02 from Wiese et al.,2019 '
da.attrs['contributions']='ANT: ice melt from Antarctica Ice Sheet (including adjacente glaciers, with 300km buffer)\
                            GRE: ice melt from Greenland Ice Sheet (inclduing djacent glaciers, with 300km buffer)\
                            GLA: ice melt from Glaciers (excluding adjacent from Greenland and Antarctica\
                            TWS: terrestrial water storage contribution (exluding glaciers and Ice Sheets)\
                            LWS: land water storage contribution (TWS+Glaciers, excluding Ice Sheets)\
                            BAR: Total baystatic contribution from land (ANT+GRE+GLA+TWS), full land signal from mascons'
da.attrs['WaterGAP_comment']='WGP_gl is an integrated model of TWS with glacier estimates from Marzeion (2012). WGP_std is the standard TWS model'
da.attrs['GIA']='Both CSR and JPL mascons had GIA signals removed by using GIA estimates rom ICE6G-D (Peltier, 2016)'
da.attrs['Mascons_comment']='Mascons had time mean from 2004-2009.999 removed'
da.attrs['months_missing']='Months Missing from Mascons: 2002-06;2002-07;2003-06;2011-01;2011-06;2012-05;2012-10;2013-03;2013-08;2013-09;2014-02;2014-07;2014-12;2015-06;2015-10;2015-11;2016-04;2016-09;2016-10;2017-02;2017-07;2017-08;2017-09;2017-10;2017-11;2017-12;2018-01;2018-02;2018-03;2018-04;2018-05;2018-08-2018-09'
da.attrs['resolution']='Data is given in 1deg resolution. Original WaterGAP, PCR-GlobWB and Mascons were in 0.5 degree. \
                        IMB, R19,M19 and ZMP were global values redistributed accordingly to each region (see mask file).'
da['mask'].attrs['long_name']='Barystatic ontributions mask'
da['mask2'].attrs['long_name']='Barystatic ontributions mask'
da['mask'].attrs['code']=ds.code_mask6
da['mask2'].attrs['code']=ds.code_mask12
da.attrs['time_coverage']='1993.01-2020.08 Datasets that do to cover the whole period have NaNs\
    R19,M19 and ZMP were originally yearly data. Other datasets have been averaged to obtain an year mean.' 
da.attrs['time_range']='1993.01-2020.08. Time is standard gregorian calendar and tdec is decimal time given in months since 1993,\
    year is in years since 1993.'
da.attrs['masked_regions']='Masked regions of Mascons have been replaced by 9999(*10)'
da.attrs['script']='barystatic_standardize.py'
da.attrs['Comment']='different sources of barystatic contributions standardized to use as input for the Sl Equation'
da.attrs['Author']='Carolina M.L. Camargo'
da.attrs['reference']='Camargo et al. 2021'
da.attrs['date_created']=str(dt.datetime.now())

da.to_netcdf(path+'comb/ALL_datasets_1993-2020_180x360_v2.nc')



