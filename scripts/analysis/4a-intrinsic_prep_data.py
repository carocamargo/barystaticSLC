#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  7 17:34:55 2022

Intrinsic Uncertainty - Prep data

@author: ccamargo
"""


import numpy as np
import xarray as xr
import sys
sys.path.append("/Users/ccamargo/Documents/py_scripts/")
import utils_SL as sl 
import utils_SLE_v2 as sle

import pandas as pd
import datetime as dt

#%%
import pickle
def load_dict(name, path):
    with open(path + name + '.pkl', 'rb') as f:
        return pickle.load(f)
    
def get_mask():
    name = 'mask'
    path = '/Volumes/LaCie_NIOZ/data/barystatic/revisions/'
    
    dic = load_dict(name,path)

    return dic
#%% open mask
dic = get_mask()
# '0:ocean, 1:TWS, 2:GLA, 3:TWS/GLA, 4:AIS, 5:GIS, 6:300km AIS filter, 7:300km Greenland filter'
regions=['ocean','tws','gla','glatws','ant','gre','antbuf','grebuf']
mask = np.array(dic['mask'])
# ds=xr.open_dataset('/Volumes/LaCie_NIOZ/data/barystatic/masks/barystatic_mask2.nc')
# ds
# # ds.gre.plot()
# # ds.mask3.plot()

# mask=np.array(ds.mask2)
# # '0=ocean,1=land,2=gre 300km buffer ,3=ant 300km buf,4=glaciers'
regions=['ocean','tws','gre','ant','gla']

#%%
path_to_save='/Volumes/LaCie_NIOZ/data/barystatic/revisions/input_data/intrinsic/'

#%% JPL mascon
print('JPL mascon - get data')
# without buffer

# maskant=np.copy(mask)
# maskant[np.where(mask==3)]=100
# maskant[np.where(maskant<100)]='nan'

#% % mascons
dataset=['JPL']
# path='/Volumes/LaCie_NIOZ/data/barystatic/original/GRA-Mascon-'
path='/Volumes/LaCie_NIOZ/data/barystatic/regrid/GRA-Mascon-'

fin=path+'JPL_180x360.nc'
ds=xr.open_dataset(fin)
# print(ds)

time=np.arange(0,len(ds.time))


for i,reg in  enumerate(regions):
    data_sl=np.zeros((len(time),len(ds.lat),len(ds.lon)))
    data_unc=np.full_like(data_sl,0)
    if i==0:
        msk = np.full_like(mask,np.nan)
        msk[mask==0] = 1
    if i==1:
        msk = np.full_like(mask,np.nan)
        msk[mask==1] = 1
        msk[mask==3] = 1
    if i==2:
        msk = np.full_like(mask,np.nan)
        msk[mask==5] = 1
        msk[mask==7] = 1
    if i==3:
        msk = np.full_like(mask,np.nan)
        msk[mask==4] = 1
        msk[mask==6] = 1
    if i==4:
        msk = np.full_like(mask,np.nan)
        msk[mask==2] = 1        
        
        
    # if d=='JPL':
    #     data_unc=np.zeros((len(time),len(ds.lat),len(ds.lon)))
        
    for j in time:
        # tmp_unc=np.array(ds.uncertainty[j,:,:]*msk)
        # tmp_sl=np.array(ds.uncertainty[j,:,:]*msk)
        
        
        data_sl[j,:,:]=np.array(ds.uncertainty[j,:,:]*msk)
        data_unc[j,:,:]=np.array(ds.uncertainty[j,:,:]*msk)
        
        # if d=='JPL':
        #     tmp_unc=np.array(ds.uncertainty[j,:,:])
        #     tmp_unc[np.where(mask!=i)]='nan'
        #     data_unc[j,:,:]=tmp_unc
    
    # out 
    if i==0:
        da=xr.Dataset(data_vars={reg+'_slc':(('time','lat','lon'),data_sl),
                                 reg+'_unc':(('time','lat','lon'),data_unc),
                                 },
                      coords={'time':ds.time,
                              'lat':ds.lat,
                              'lon':ds.lon})
        
    da[reg+'_slc']=(('time','lat','lon'),data_sl)
    da[reg+'_unc']=(('time','lat','lon'),data_unc)
    
    # if d=='JPL':
    #     ds[reg+'_unc']=(('time','lat','lon'),data_unc)    
    
da['mask']=(('lat','lon'),mask)
da.attrs['regions']='Regions selected from the mascons using the barystatic component mask'
da.attrs['bary_mask']='0=ocean,1=land,2=gre 300km buffer ,3=ant 300km buf,4=glaciers'
da.attrs['script']='4a_intrinsic_prep_data.py'
da.attrs['missing_value']='masked areas were replaced by nans'

#%% make dates
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
    
#%% standardize time
names=regions[1:len(regions)]
# ds=xr.open_dataset(f)
da_y=da.groupby('time.year').mean(dim='time')

#% % make empty arrays:
SL_mm=np.zeros((len(names),len(time),180,360))
SL_mm.fill('nan')
SL_EWH=np.full_like(SL_mm,np.nan)
UNC_EWH =np.full_like(SL_mm,np.nan)
UNC_mm=np.full_like(SL_mm,np.nan) 

SL_mm_y=np.zeros((len(names),len(timey),180,360))
SL_mm_y.fill('nan')
SL_EWH_y=np.full_like(SL_mm_y,np.nan)
UNC_mm_y=np.full_like(SL_mm_y,np.nan)
UNC_EWH_y=np.full_like(SL_mm_y,np.nan)

# % % find time index
time_local=np.array(ds.time)
tdec_local,tdec0=sl.get_dec_time(time_local,ns=1e-9)
lon=np.array(da.lon)
lat=np.array(da.lat)
timey_local=np.array(da_y.year)

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

# #############
# loop over variables 
for ivar, var in enumerate(names):
    # # Water thickness
    data_sl=np.array(da[var+'_slc'])*10 # mm of water thickness
    datay_sl=np.array(da_y[var+'_slc'])*10 
    data_unc=np.array(da[var+'_unc'])*10 # mm of water thickness
    datay_unc=np.array(da_y[var+'_unc'])*10 
    
    # Put data in main array
    SL_EWH_y[ivar,idx_0y:idx_1y,:,:]=datay_sl[idx_00y:len(timey_local),:,:]
    UNC_EWH_y[ivar,idx_0y:idx_1y,:,:]=datay_unc[idx_00y:len(timey_local),:,:]
    
    for i,j in enumerate(idx_local_m):
        if np.isfinite(j):
            j=int(j)
            SL_EWH[ivar,i,:,:]=data_sl[j,:,:]
            UNC_EWH[ivar,i,:,:]=data_unc[j,:,:]
            
    # Transform data: 
    for i in range(len(data_sl)):
        data_sl[i,:,:]=sle.EWH_to_height(data_sl[i,:,:]).reshape(180,360)
        data_unc[i,:,:]=sle.EWH_to_height(data_unc[i,:,:]).reshape(180,360)
        
    for i in range(len(datay_sl)):
        datay_sl[i,:,:]=sle.EWH_to_height(datay_sl[i,:,:]).reshape(180,360)
        datay_unc[i,:,:]=sle.EWH_to_height(datay_unc[i,:,:]).reshape(180,360)
        
    # Put transformed data in main array:
    SL_mm_y[ivar,idx_0y:idx_1y,:,:]=datay_sl[idx_00y:len(timey_local),:,:]
    UNC_mm_y[ivar,idx_0y:idx_1y,:,:]=datay_unc[idx_00y:len(timey_local),:,:]
    for i,j in enumerate(idx_local_m):
        if np.isfinite(j):
            j=int(j)
            SL_mm[ivar,i,:,:]=data_sl[j,:,:]
            UNC_mm[ivar,i,:,:]=data_unc[j,:,:]
   
#%% make and save data array

#% %
da2=xr.Dataset(data_vars={'SL_mm':(('name','time','lat','lon'),SL_mm),
                          'SL_EWH':(('name','time','lat','lon'),SL_EWH),
                          'SL_mm_y':(('name','year','lat','lon'),SL_mm_y),
                          'SL_EWH_y':(('name','year','lat','lon'),SL_EWH_y),
                          'UNC_mm':(('name','time','lat','lon'),UNC_mm),
                          'UNC_EWH':(('name','time','lat','lon'),UNC_EWH),
                          'UNC_mm_y':(('name','year','lat','lon'),UNC_mm_y),
                          'UNC_EWH_y':(('name','year','lat','lon'),UNC_EWH_y),
                          
                          # 'mask':(('lat','lon'),ds_mask['mask6']),
                          # 'mask2':(('lat','lon'),ds_mask['mask12']),
                         
                        },
                          coords={'lat':ds.lat,
                                  'lon':ds.lon,
                                  'time':time2,
                                  'tdec':tdec,
                                  'year':timey,
                                  'name':names})
da2['SL_mm'].attrs['units']='mm of sea level'
da2['SL_mm_y'].attrs['units']='mm of sea level'
da2['SL_EWH'].attrs['units']='mm of Equivalent Water Thickness'
da2['SL_EWH_y'].attrs['units']='mm of Equivalent Water Thickness'
da2['SL_mm'].attrs['long_name']='Monthly ocean mass change in mm of sea level height'
da2['SL_EWH'].attrs['long_name']='Monthly ocean mass change in mm of equivalent water thickness'
da2['SL_mm_y'].attrs['long_name']='Yearly averages of ocean mass change in mm of sea level height'
da2['SL_EWH_y'].attrs['long_name']='Yearly averages of ocean mass change in mm of equivalent water thickness'
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
da2['name'].attrs['long_name']='Name of the source region/contribution of ocean mass change'
#% %
da2.attrs['metadata']='Ocean mass sea level changes (in mm and EWH), from differente contributions (ANT, GRE, TWS, GLA) \
      from JPL Mascons)'

da2.attrs['Author']='Carolina M.L. Camargo'
da2.attrs['reference']='Camargo et al. 2021'
da2.attrs['date_created']=str(dt.datetime.now())

da2.to_netcdf(path_to_save+'JPL.nc')


djpl = xr.open_dataset(path_to_save+'JPL.nc')
dimtime=len(djpl.tdec)
#%% IMBIE
dataset='GIS_IMB'
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

df = df.set_index('Year')
df = df.dropna() # From 1992 onwards

gre_regions = ['GIS','SMB','DYN']
cols=[2,6,10]
slc=np.zeros((len(gre_regions),len(df)))
unc=np.full_like(slc,0)

for i,col in enumerate(cols):
    # print(df.columns[col])
    slc[i,:]=np.array(df[df.columns[col]])
    unc[i,:]=np.array(df[df.columns[col+1]])

da=xr.Dataset(data_vars={'{}_slc'.format(dataset):(('reg_{}'.format(dataset),'time'),slc),
                        '{}_unc'.format(dataset):(('reg_{}'.format(dataset),'time'),unc),                           
                        },
                          coords={
                                  'time':np.array(df.index),
                                  'reg_{}'.format(dataset):gre_regions
                                  })

time_dim=len(df)
#% % get ANT dataset
# IMBIE 2018 Antartic Dataset
#This spreadsheet contains the IMBIE-2018 datasets for Antarctica (worksheet 1),
#Antarctic Peninsula (worksheet 2), East Antarctica (worksheet 3) and West Antarctica (worksheet 4). 
#Each worksheet includes data on monthly cumulative ice sheet mass changes and their estimated uncertainty. 
#The data are expressed in units of mass (Gigatons – columns B and C) 
#and in units of equivalent mean global sea level rise (millimetres – columns D and E).

dataset='AIS_IMB'
file='/Volumes/LaCie_NIOZ/data/barystatic/original/ANT_IMBIE.xlsx'
regions = ['AIS','AP','EAIS','WAIS']
local_dim = len(pd.read_excel(file))
unc=np.zeros((len(regions),time_dim))
slc = np.full_like(unc,0)

for sheet, reg in enumerate(regions):
    df=pd.read_excel(file,
                  #sheet_name=[0,1] # load first and second sheet as a dict of df
                  sheet_name=sheet # sheet 1
                  )
    unc[sheet,0:local_dim] = np.array(df['Cumulative sea level contribution uncertainty (mm)'])
    slc[sheet,0:local_dim] = np.array(df['Cumulative sea level contribution (mm)'])

da['{}_slc'.format(dataset)]=(('reg_{}'.format(dataset),'time'),slc)
da['{}_unc'.format(dataset)]=(('reg_{}'.format(dataset),'time'),unc)

da = da.assign_coords({'reg_{}'.format(dataset):regions})

#%%
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
#%% open mask
name_mask = 'masks_dict'
path_mask = '/Volumes/LaCie_NIOZ/data/barystatic/'
# path_mask=path+name
masks = load_dict(name_mask,path_mask)
print(masks.keys())
mask_keys = {'AIS_IMB':'AIS_regions',
                 'AIS_UCI':'AIS_basins',
                 'GIS_IMB':'GIS_regions',
                 'GIS_UCI':'GIS_basins',
                 'GLA_ZMP':'Glaciers'
                 }
#%% 
dimlat=180;dimlon=360
# dimtime=len(da.time)
names=['GIS_IMB','AIS_IMB']
for iname,name in enumerate(names):
    masks = load_dict(name_mask,path_mask)
    regional_slc_mm = np.zeros((len(da['reg_'+name]),dimtime,dimlat,dimlon))
    regional_slc_mm.fill('nan')
    regional_unc_mm = np.full_like(regional_slc_mm,np.nan)
    regional_slc_ewh = np.full_like(regional_slc_mm,np.nan)
    regional_unc_ewh = np.full_like(regional_slc_mm,np.nan)

            
    for ireg,reg in enumerate(np.array(da['reg_'+name])):
        print(reg)
#% %
        mask = masks[mask_keys[name]][reg]
        for itime,t in enumerate(da.time):
            regional_slc_mm[ireg,itime,:,:] = glb_to_reg(
                                    np.array(da[name+'_slc'][ireg,itime]),mask)
            regional_unc_mm[ireg,itime,:,:] = glb_to_reg(
                                    np.array(da[name+'_unc'][ireg,itime]),mask)
            regional_slc_ewh[ireg,itime,:,:] = sle.height_to_EWH(glb_to_reg(
                                    np.array(da[name+'_slc'][ireg,itime]),mask)).reshape(180,360)
            regional_unc_ewh[ireg,itime,:,:] = sle.height_to_EWH(glb_to_reg(
                                    np.array(da[name+'_unc'][ireg,itime]),mask)).reshape(180,360)
        
    # make regional dataset:
    da2=xr.Dataset(data_vars={
        'SL_mm':(('reg','time','lat','lon'),regional_slc_mm),
        'UNC_mm':(('reg','time','lat','lon'),regional_unc_mm),
        'SL_EWH':(('reg','time','lat','lon'),regional_slc_ewh),
        'UNC_EWH':(('reg','time','lat','lon'),regional_unc_ewh),
        },
                      coords={'lon':np.arange(0.5,360,1),
                              'lat':np.arange(-89.5,90,1),
                              'reg':[str(reg) for reg in np.array(da['reg_'+name])],
                              'time':djpl.time
                              })
    
    reg=name.split('_')[0]
    da2=da2.drop_sel(reg=reg)
    mask=masks[reg+'_regions'][reg]
    da2=da2.sum(dim='reg')* mask
    # save 
        
    da2.to_netcdf(path_to_save+name+'.nc')
            



#%% open datasets and merge IMBIE and JPL
print('combine IMBIE and JPL')
names=['LWS_JPL','GIS_JPL','AIS_JPL','GLA_JPL',
       'AIS_IMB','GIS_IMB']
flist=[path_to_save+'JPL.nc',
       path_to_save+'AIS_IMB.nc',
       path_to_save+'GIS_IMB.nc',
       ]

 
#% % combine

SL_mm=np.zeros((len(names),dimtime,180,360))
SL_mm.fill('nan')
SL_EWH=np.full_like(SL_mm,np.nan)
UNC_EWH=np.full_like(SL_mm,np.nan)
UNC_mm=np.full_like(SL_mm,np.nan)


for ifile, f in enumerate(flist):
    ds=xr.open_dataset(f)
    if ifile==0:
        j=len(ds.name)
        SL_mm[ifile:j,:,:,:]=np.array(ds.SL_mm)
        SL_EWH[ifile:j,:,:,:]=np.array(ds.SL_EWH)
        UNC_mm[ifile:j,:,:,:]=np.array(ds.UNC_mm)
        UNC_EWH[ifile:j,:,:,:]=np.array(ds.UNC_EWH)
        j=j-1
    else:
        ifile=j+1
        j=j+1
        # print(ifile)
        # j=len(ds.name)+j
        SL_mm[ifile,:,:,:]=np.array(ds.SL_mm)
        SL_EWH[ifile,:,:,:]=np.array(ds.SL_EWH)
        UNC_mm[ifile,:,:,:]=np.array(ds.UNC_mm)
        UNC_EWH[ifile,:,:,:]=np.array(ds.UNC_EWH)

#% % save 
da=xr.Dataset(data_vars={'SL_mm':(('name','time','lat','lon'),SL_mm),
                           'SL_EWH':(('name','time','lat','lon'),SL_EWH),
                           'UNC_mm':(('name','time','lat','lon'),UNC_mm),
                           'UNC_EWH':(('name','time','lat','lon'),UNC_EWH),

                          
                         },
                           coords={'lat':ds.lat,
                                   'lon':ds.lon,
                                   'time':ds.time,
                                   'tdec':djpl.tdec,
                                   'name':names})

da['SL_mm'].attrs['units']='mm of sea level'
# da['SL_mm_y'].attrs['units']='mm of sea level'
da['SL_EWH'].attrs['units']='mm of Equivalent Water Thickness'
# da['SL_EWH_y'].attrs['units']='mm of Equivalent Water Thickness'
da['SL_mm'].attrs['long_name']='Monthly ocean mass change in mm of sea level height'
da['SL_EWH'].attrs['long_name']='Monthly ocean mass change in mm of equivalent water thickness'
# da['SL_mm_y'].attrs['long_name']='Yearly averages of ocean mass change in mm of sea level height'
# da['SL_EWH_y'].attrs['long_name']='Yearly averages of ocean mass change in mm of equivalent water thickness'

da.attrs['script']='intrinsic_unc-SLE.py'
da.attrs['Comment']='different sources of barystatic contributions standardized to use as input for the Sl Equation'
da.attrs['Author']='Carolina M.L. Camargo'
da.attrs['reference']='Camargo et al. 2021'
da.attrs['date_created']=str(dt.datetime.now())

# path='/Volumes/LaCie_NIOZ/data/barystatic/intrinsic_unc/use/'+'comb/'
da.to_netcdf(path_to_save+'SL_unc_JPL_IMBIE.nc')
