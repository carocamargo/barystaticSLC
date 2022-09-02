#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  7 17:48:16 2022
Compute intrinsic unc rate and fingeprint
Update: Monte-Carlo perturbation of the trend to 
obtain intrinsic unc

@author: ccamargo
"""


import numpy as np
import xarray as xr
import sys
sys.path.append("/Users/ccamargo/Documents/py_scripts/")
import utils_SL as sl 
import utils_SLE_v2 as sle
sys.path.append("/Users/ccamargo/Documents/py_scripts/OM/")
import utils_OM as om 
from datetime import datetime as dt

#%%
path='/Volumes/LaCie_NIOZ/data/barystatic/revisions/input_data/intrinsic/'
# path='/Volumes/LaCie_NIOZ/data/barystatic/intrinsic_unc/use/'+'comb/'
infile='SL_unc_JPL_IMBIE.nc'
# print(ds)
# print(ds.name)

# select time
# tdec=np.array(ds.tdec)
periods =[ 
        # (2005,2016),
        (1993,2018),
           # (1993,2017),
        (2003,2017)
          ]
period=periods[-1]
#% %
print('Computing intrinsic trends\n')

for period in periods:
    start = dt.now()
    #% %
    print('period: {}'.format(period))
    # period=periods[-1]
    ds=xr.open_dataset(path+infile)
    # ds_sl=xr.open_dataset('/Volumes/LaCie_NIOZ/data/barystatic/use/comb/ALL_datasets_1993-2020_180x360_v3_update.nc')
    t0=period[0]
    t1=period[1]-1
    ds['tyear']=(('time'),ds.tdec)
    # Select time
    to=str(t0)+'-01-01';ti=str(t1)+'-01-01'
    ds= ds.sel(time=slice(to,ti))
    # ds_sl= ds_sl.sel(time=slice(to,ti))
    
    tdec=np.array(ds.tyear)
    lat=np.array(ds.lat)
    lon=np.array(ds.lon)
    names=np.array(ds.name)
    # names = [name for name in names if not name=='TCWS_JPL']
    if t0< 2002: # before grace
        names = [name for name in names if not name.endswith('JPL')]
    # names_sl = [name for name in np.array(ds_sl.name) ]
    # for i,name in enumerate(names_sl):
    #     if '300' in name: 
    #         name=name.split('_')[0]+'_'+name.split('_')[-1]
    # for name in names: 
    #     if name  in names_sl: 
    #         print(name)
    # ds_sl['name']=names_sl
    intr_unc= np.zeros((len(names),len(ds.lat),len(ds.lon)))
    trend = np.zeros((len(names),len(ds.lat),len(ds.lon)))
    trend2 = np.zeros((len(names),len(ds.lat),len(ds.lon)))
    std_trend = np.zeros((len(names),len(ds.lat),len(ds.lon)))
    
    
    slf=np.zeros((len(names),len(ds.lat),len(ds.lon)))
    slf2=np.zeros((len(names),len(ds.lat),len(ds.lon)))
    
    
    for iname,name in enumerate(names):
        #% %
        # iname=0;name=names[iname]
        da=ds.sel(name=name)
        # da_sl=ds_sl.sel(name=name)
        
        print(name)
        # select dataset
        unc=np.array(da.UNC_mm[:,:,:])
        y=np.array(da.SL_mm[:,:,:])
        # y_up = np.array(y+unc)
        # y_low = np.array(y-unc)
        # compute rate of intrinsic uncertainty
        intr_unc[iname,:,:] = om.intrinsic(time=tdec,rate=y, unc = unc, 
                                           N=1000,ci_level=0.95)
        print('{} /dataset'.format(dt.now()-start))
        trend[iname,:,:], _,_,trend2[iname,:,:],std_trend[iname,:,:]=sl.get_OLS_trend(tdec, y, sigma=unc,lat=lat, lon=lon)
        # trend_up, _=sl.get_reg_trend_OLS(tdec, y_up, lat, lon)
        # trend_low, _=sl.get_reg_trend_OLS(tdec, y_low, lat, lon)
        # trend_bound[iname,:,:] = np.max([trend_up-trend, trend-trend_low],axis=0)

        slf2[iname,:,:] = sle.run_SLE(sle.height_to_EWH(std_trend[iname,:,:]).reshape(180,360),name)
        slf[iname,:,:] = sle.run_SLE(sle.height_to_EWH(intr_unc[iname,:,:]).reshape(180,360),name)
        
        # ensure that the uncertainties are combined in quadrature for the fingerprint
        # slf[iname,:,:] = sle.run_SLE(sle.height_to_EWH(intr_unc[iname,:,:]*intr_unc[iname,:,:]).reshape(180,360),name)
        # slf[iname,:,:]=np.sqrt(np.abs(slf[iname,:,:]))
        
        # slf2[iname,:,:] = sle.run_SLE(sle.height_to_EWH(std_trend[iname,:,:]*std_trend[iname,:,:]).reshape(180,360),name)
        # slf2[iname,:,:]=np.sqrt(np.abs(slf2[iname,:,:]))
        
    #% % make data array
    da=xr.Dataset(data_vars={'intrinsic_unc_source':(('name','lat','lon'),intr_unc), 
                             'intrinsic_unc_SLF':(('name','lat','lon'),slf), 
                             'intrinsic_unc_source_old':(('name','lat','lon'),std_trend), 
                             'intrinsic_unc_SLF_old':(('name','lat','lon'),slf2), 
                                                      
    
                             },
                               coords={'lat':lat,
                                       'lon':lon,
    
                                       'name':names})
    da.attrs['units']='mm/year'
    da['lat'].attrs['standard_name']='latitude'
    da['lat'].attrs['long_name']='Latitude'
    da['lat'].attrs['units']='degrees_north'
    da['lat'].attrs['axis']='Y'
    
    da['lon'].attrs['standard_name']='longitude'
    da['lon'].attrs['long_name']='Longitude'
    da['lon'].attrs['units']='degrees_east'
    da['lon'].attrs['axis']='X'
    
    da.attrs['metadata']='Intrinsic uncertainty obtained by computing the standard deviation of the trend, when propagating the unc in the OLS'
    da.attrs['method']="perturb each estimate 1000 times, and take the 95% CI of the estimtes distribution"



    # % %  save 
    pwd='/Volumes/LaCie_NIOZ/data/barystatic/revisions/results_final/'
    da.to_netcdf(pwd+'{}-{}/intrinsic_unc_prop_{}-{}.nc'.format(t0,t1,t0,t1))
    # da.to_netcdf(path+'source_prop_SLF_{}-{}.nc'.format(t0,t1))
    
# %% plot
# from cartopy import crs as ccrs , feature as cfeature
# import matplotlib.pyplot as plt
# t0=2003
# t1=2016
# pwd='/Volumes/LaCie_NIOZ/data/barystatic/revisions/results/'
# da= xr.open_dataset(pwd+'{}-{}/intrinsic_unc_prop_{}-{}.nc'.format(t0,t1,t0,t1))
# names = np.array(da.name)
# landcolor='darkgrey' 
# dpi=300
# clim=0.5
# cmap='Blues'
# cmin=0
# cmax=clim
# interval = 0.1
# X=np.array(da.lon)
# Y=np.array(-da.lat)
# fontsize=25
# ticksize=20
# fig = plt.figure(figsize=(15,10), facecolor='w',dpi=dpi)

# for iname,name in enumerate(names):
#     da2=da.sel(name=name)
#     # print(name)
#     ax1 = plt.subplot(3,2,iname+1, projection=ccrs.Robinson())

#     # SLF monthly then trend:
#     data=np.array(np.abs(da2.intrinsic_unc_SLF [:,:]))
# #    print(np.nanmin(data))
# #    print(np.nanmax(data))
#     # ax1=plt.subplot(111,projection=ccrs.Robinson())
#     ax1.coastlines(resolution='110m', zorder=3,color=landcolor) # zorder=3 makes sure that no other plots overlay the coastlines
#     ax1.add_feature(cfeature.LAND,color=landcolor,# alpha=0.5,
#                       zorder=3)
#     ax1.set_global() # make sure the projection is maximised inside the plot. Makes a circle

#     # make it discrete
  
#     lv=np.arange(cmin,cmax+interval,interval)
#     csf=plt.contourf(X,Y,np.abs(data),levels=lv,
#               transform = ccrs.PlateCarree(),cmap=cmap)

#     mu,glb,mask=sle.reg_to_glb(data,Y,X)
#     print(mu)
#     glb=np.round(glb,3)
 
#     cs=ax1.contour(X,Y,data,
#                     levels=[glb],
#                 # levels=[0.15, 0.30, 0.45, 0.6, 0.75],
#             # vmin=-0.6,
#             # vmax=0.6,
#         transform = ccrs.PlateCarree(),
#         #cmap='coolwarm',#extend='both'
#         colors=('black',),linestyles=('--',),linewidths=(2,)
#         )   
#     ax1.clabel(cs,cs.levels,fmt='%5.2f',colors='k',fontsize=12)

#     cp=plt.pcolormesh(X,Y,np.abs(data),
#                 vmin=cmin,vmax=cmax,
#                 zorder=0,
#                 transform = ccrs.PlateCarree(),cmap=cmap)
#     ax1.set_title('({}). '.format(str(name)),size=fontsize)   


# cbar_ax2 = fig.add_axes([0.152, 0.05, 0.72, 0.033])
# cbar2 = plt.colorbar(csf, cax=cbar_ax2, orientation='horizontal')
# cbar2.set_label(label='Intrinsic Uncertainty \n{}-{} (mm/yr)'.format(t0,t1),size=fontsize, family='serif')
# cbar2.ax.tick_params(labelsize=ticksize) 

# plt.show()
    
