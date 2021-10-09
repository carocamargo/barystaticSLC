#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 28 08:25:15 2021

Noise model selection 

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
#%%

periods =[ (2005,2016),(1993,2018),# 
          (2003,2017)
          ]
for period in periods:
    #% %
    t0=period[0]
    t1=period[1]-1
    
    # pwd = '/Volumes/LaCie_NIOZ/data/barystatic/hector/source/local_test/'
    # pwd = '/Volumes/LaCie_NIOZ/data/barystatic/hector/source/EWH/'
    
    pwd = '/Volumes/LaCie_NIOZ/data/barystatic/hector/source/mixed/'
    
    # t0=2005
    # t1=2015
    pwd = pwd+'{}-{}/'.format(t0,t1)
    
    flist=sl.get_filelist(pwd,'*.nc')
    
    flist=[file for file in flist if not file.split('_')[4]=='unc']
    flist=[file for file in flist if not file.split('_')[3]=='unc']
    
    ds=xr.open_dataset(flist[2])
    print(ds)
    for file in flist:
        ds=xr.open_dataset(file)
        
        region = file.split('/')[-1].split('_')[0]
        dataset = file.split('/')[-1].split('_')[1]
        if dataset == '300':
            dataset = file.split('/')[-1].split('_')[2]
        name = '{}_{}'.format(region,dataset)
        
        if len(ds.trend.shape) ==3: # regional data (nm, lat, lon)
            #% % Rank noise models for each dataset
            lat=np.array(ds.lat)
            lon=np.array(ds.lon)
            trend=np.array(ds.trend)
            uncert=np.array(ds.unc)
            mask=np.array(trend[0,:,:])
            mask[np.isfinite(mask)]=1
            aic=np.array(ds.aic)
            bic=np.array(ds.bic)
            bic_tp=np.array(ds.bic_tp)#np.nanmean([aic,bic],axis=0)
            nm=np.array(ds.nm)
            
            #% %
            rkb=np.full_like(np.zeros((len(ds.lat),len(ds.lon))),np.nan)
            rka=np.full_like(rkb,np.nan)
            rk=np.full_like(rkb,np.nan)
            
            best_trend_aic=np.full_like(rkb,np.nan)
            best_trend_bic=np.full_like(rkb,np.nan)
            best_trend_bic_tp=np.full_like(rkb,np.nan)
            best_unc_aic=np.full_like(rkb,np.nan)
            best_unc_bic=np.full_like(rkb,np.nan)
            best_unc_bic_tp=np.full_like(rkb,np.nan)
            
            
            
            ilat= np.where(np.isfinite(mask))[0][0]
            ilon=  np.where(np.isfinite(mask))[1][0]
            for ilat,l in enumerate(lat):
                for ilon,l in enumerate(lon):
                    
                    if np.isfinite(mask[ilat,ilon]) and np.any(np.isfinite(aic[:,ilat,ilon])):
                        # # Quick QC
                        # std = np.nanstd(trend[ilat,ilon]) + np.nanmax(uncert[ilat,ilon])
                        # mu = np.nanmean(trend[ilat,ilon])
                        # if np.any(trend[ilat,ilon]> mu +std) or np.any(trend[ilat,ilon]< mu-std):
                        #     idx=np.array(np.where((trend[ilat,ilon] > mu+std) | (trend[ilat,ilon] < mu-std)))[0]
                        #     for ix in idx:
                        #         # print('QC removed {}'.format(nm[ix]))
                        #         trend[jname,ix,ilat,ilon]=np.nan;
                        #         uncert[jname,ix,ilat,ilon]=np.nan
                        #         # aic[jname,ix,ilat,ilon]=np.nan;
                        #         # bic[jname,ix,ilat,ilon]=np.nan;
                        #         # logL[jname,ix,ilat,ilon]=np.nan
                        #         aic[jname,ix,ilat,ilon]=np.nanmax(aic[ilat,ilon])
                        #         bic[jname,ix,ilat,ilon]=np.nanmax(bic[ilat,ilon])
                        #         bic_tp[jname,ix,ilat,ilon]=np.nanmax(bic_tp[ilat,ilon])
                        #         logL[jname,ix,ilat,ilon]=np.nanmax(logL[ilat,ilon])
                                
                        targeta=aic[:,ilat,ilon] # All datasets and noise-models for this lat and lon
                        targetb=bic[:,ilat,ilon]       
                        targetab=bic_tp[:,ilat,ilon]
                        
                        bic2=np.zeros((len(nm)))
                        aic2=np.zeros((len(nm)))
                        bic_tp2=np.zeros((len(nm)))
                
                        # for iname in np.arange(0,len([name])):
                        iname = 0
            ##                AIC
                        tg=targeta[:] # all noise models for this dataset
                        logaic=np.zeros((tg.shape))
                        for inm in np.arange(0,len(nm)):
                            logaic[inm]=np.exp((np.nanmin(tg)-tg[inm])/2)
                        
                        logaic 
                        indaic=np.where(logaic>0.5)
                        indaic=np.array(indaic)
                        for inm in np.arange(0,len(nm)):
                            for ind in indaic:
                                for ind2 in ind:
                                    if inm==ind2:        
                                        aic2[inm]=1
            #                
                        
                        #BIC
                        tg=targetb[:]
                        tg
                        logbic=np.zeros((tg.shape))
                        for inm in np.arange(0,len(nm)):
                            logbic[inm]=np.exp((np.nanmin(tg)-tg[inm])/2)
                            # logbic[inm] = tg[inm] - np.nanmin(tg)
                        
                        logbic 
                        indbic=np.where(logbic>0.5)
                        # indbic=np.where(logbic<2)
                        indbic=np.array(indbic)
                        indbic
                        for inm in np.arange(0,len(nm)):
                            for ind in indbic:
                                for ind2 in ind:
                                    if inm==ind2:        
                                        bic2[inm]=1
            #                AIC|BIC
                        tg=targetab[:]
                        tg
                        logab=np.zeros((tg.shape))
                        for inm in np.arange(0,len(nm)):
                            logab[inm]=np.exp((np.nanmin(tg)-tg[inm])/2)
                        
                        logab 
                        indab=np.where(logab>0.5)
                        indab=np.array(indab)
                        indab
                        for inm in np.arange(0,len(nm)):
                            for ind in indab:
                                for ind2 in ind:
                                    if inm==ind2:        
                                        bic_tp2[inm]=1
                                        
                    
                        # indaic=np.sum(aic2,axis=0)                   
                        # tmp=np.array(np.where(aic2==np.nanmax(indaic)))
                        inda= np.where(aic2==np.max(aic2))
                        tmp=np.array(np.where(targeta==targeta[inda].min()))
                        rka[ilat,ilon]=tmp[0,0]
            
                        # indbic=np.sum(bic2,axis=0)
                        # tmp=np.array(np.where(bic2==np.nanmax(indbic)))
                        indb= np.where(bic2==np.max(bic2))
                        tmp=np.array(np.where(targetb==targetb[indb].min()))
                        rkb[ilat,ilon]=tmp[0,0]
            
                        # indab=np.sum(bic_tp2,axis=0)
                        indab= np.where(bic_tp2==np.max(bic_tp2))
                        tmp=np.array(np.where(targetab==targetab[indab].min()))
                        rk[ilat,ilon]=tmp[0,0]
            
                        
            #% %
            rka05=np.array(rka)
            rkb05=np.array(rkb)
            rk05=np.array(rk)
            trend05=np.array(trend)
            unc05=np.array(uncert)
            
            # Compute % of grid cells and % ocean area for each noise model, 
            # according to each selection criterea and time period
            
            data=rka05[:,:]
            
            perc=np.zeros((len(nm)))
            # count nans:
            tmp=list(np.hstack(data))
            j=np.array(np.where(np.isnan(tmp)))
            
            
            for i in range(len(nm)):
               # print(i)
                #print(NM[i])
                tmp=list(np.hstack(data))
                perc[i]=(tmp.count(i)/(len(tmp)-len(j[0,:])))*100
                #print(perc_bic[i])
            
            rel_area=np.zeros((len(nm)))
            dry_msk=sl.get_dry_msk(data)
            total_area,area=sl.get_ocean_area(lat,lon,dry_msk,info=False)
            for i in range(len(nm)):
                tmp=np.copy(data)
                tmp[np.where(tmp!=i)]='nan'
                tmp[np.where(np.isnan(tmp)==False)]=1
                tmp[np.where(np.isnan(tmp))]=0
                surf,area=sl.get_ocean_area(lat,lon,tmp,info=False)
                rel_area[i]=(surf/total_area)*100
            
            df=pd.DataFrame(rel_area,columns=['aic_area05'])
            df['aic_grid05']=perc
            
            #% %
            data=rkb05[:,:]
            perc=np.zeros((len(nm)))
            # count nans:
            tmp=list(np.hstack(data))
            j=np.array(np.where(np.isnan(tmp)))
            
            
            for i in range(len(nm)):
               # print(i)
                #print(NM[i])
                tmp=list(np.hstack(data))
                perc[i]=(tmp.count(i)/(len(tmp)-len(j[0,:])))*100
                #print(perc_bic[i])
            
            rel_area=np.zeros((len(nm)))
            dry_msk=sl.get_dry_msk(data)
            total_area,area=sl.get_ocean_area(lat,lon,dry_msk,info=False)
            for i in range(len(nm)):
                tmp=np.copy(data)
                tmp[np.where(tmp!=i)]='nan'
                tmp[np.where(np.isnan(tmp)==False)]=1
                tmp[np.where(np.isnan(tmp))]=0
                surf,area=sl.get_ocean_area(lat,lon,tmp,info=False)
                rel_area[i]=(surf/total_area)*100
            
            df['bic_tparea05']=rel_area
            df['bic_tpgrid05']=perc
            
            #% %
            data=rk05[:,:]
            perc=np.zeros((len(nm)))
            # count nans:
            tmp=list(np.hstack(data))
            j=np.array(np.where(np.isnan(tmp)))
            
            
            for i in range(len(nm)):
               # print(i)
                #print(NM[i])
                tmp=list(np.hstack(data))
                perc[i]=(tmp.count(i)/(len(tmp)-len(j[0,:])))*100
                #print(perc_bic[i])
            
            rel_area=np.zeros((len(nm)))
            dry_msk=sl.get_dry_msk(data)
            total_area,area=sl.get_ocean_area(lat,lon,dry_msk,info=False)
            for i in range(len(nm)):
                tmp=np.copy(data)
                tmp[np.where(tmp!=i)]='nan'
                tmp[np.where(np.isnan(tmp)==False)]=1
                tmp[np.where(np.isnan(tmp))]=0
                surf,area=sl.get_ocean_area(lat,lon,tmp,info=False)
                rel_area[i]=(surf/total_area)*100
            
            df['bic_tp_area05']=rel_area
            df['bic_tp_grid05']=perc
            
            
            #% % plot preferred Noise model
            fig=plt.figure(figsize=(20,10),dpi=300)
            gs = GridSpec(3, 3, figure=fig)
            
            cmax=7
            cmin=0
            cmap=plt.cm.get_cmap('Set3', cmax+1)
            
            # AIC
            ax1 = plt.subplot(gs.new_subplotspec((0, 0),colspan=3))
            data=rka05[:,:]
            
            perc=np.array(df.aic_area05)
            
            lb=list(nm)
            for i in range(len(nm)):
                lb[i]=str(nm[i]+":"+str(np.round(perc[i],1)) )
            
            
            
            m=Basemap(projection='cyl',
                      llcrnrlon=np.min(lon),urcrnrlon=np.max(lon), 
                        llcrnrlat=lat.min(),urcrnrlat=lat.max(), 
                        resolution='c')       
            cax=m.pcolormesh(lon,lat,data,shading='flat',cmap=cmapnm)
            plt.clim(cmin,cmax)
            m.drawparallels([0])
            m.drawmeridians([180])
            
            deg=np.arange(0,360,60)
            lblon=[None]*len(deg)
            for i,d in enumerate(deg): lblon[i]=lon2str(d)
            ax1.set_xticks(deg)
            ax1.set_xticklabels(lblon,fontsize=15)
            
            deg=np.arange(-60,61,30)
            lblat=[None]*len(deg)
            for i,d in enumerate(deg): lblat[i]=lat2str(d)
            ax1.set_yticks(deg)
            ax1.set_yticklabels(lblat,fontsize=15)
            
            cbar=plt.colorbar(cax,ticks=range(len(nm)), #shrink=0.7,
                              #label='Preferred Noise Model',fontsize=15,
                              orientation='vertical')
            cbar.ax.set_yticklabels(lb,fontsize=15) # horizontally oriented colorbar
            plt.title('{} \n{}-{}'.format(name,t0,t1),fontsize=18)
            
            ax1.text(0.08,0.90,' a ',weight='bold',fontsize=18,
                    horizontalalignment='left',
                    transform=ax1.transAxes)
            plt.ylabel('AIC\n',fontsize=18)
            
            
            #%
            #########################
            # BIC
            ax1 = plt.subplot(gs.new_subplotspec((1, 0),colspan=3))
            
            data=rkb05[:,:]
            perc=np.array(df.bic_tparea05)
            lb=list(nm)
            for i in range(len(nm)):
                lb[i]=str(nm[i]+":"+str(np.round(perc[i],1)) )
                
            
            m=Basemap(projection='cyl',
                      llcrnrlon=np.min(lon),urcrnrlon=np.max(lon), 
                        llcrnrlat=lat.min(),urcrnrlat=lat.max(), 
                        resolution='c')       
            cax=m.pcolormesh(lon,lat,data,shading='flat',cmap=cmapnm)
            plt.clim(cmin,cmax)
            m.drawparallels([0])
            m.drawmeridians([180])
            
            deg=np.arange(0,360,60)
            lblon=[None]*len(deg)
            for i,d in enumerate(deg): lblon[i]=lon2str(d)
            ax1.set_xticks(deg)
            ax1.set_xticklabels(lblon,fontsize=15)
            
            deg=np.arange(-60,61,30)
            lblat=[None]*len(deg)
            for i,d in enumerate(deg): lblat[i]=lat2str(d)
            ax1.set_yticks(deg)
            ax1.set_yticklabels(lblat,fontsize=15)
            
            cbar=plt.colorbar(cax,ticks=range(len(nm)), #shrink=0.7,
                              #label='Preferred Noise Model',fontsize=15,
                              orientation='vertical')
            cbar.ax.set_yticklabels(lb,fontsize=15) # horizontally oriented colorbar
            plt.ylabel('BIC\n',fontsize=18)
            ax1.text(0.08,0.90,' c ',weight='bold',fontsize=18,
                    horizontalalignment='left',
                    transform=ax1.transAxes)
            #plt.ylabel('Latitude\n',fontsize=15)
            
            
            ####################
            # AIC|BIC
            ax1 = plt.subplot(gs.new_subplotspec((2, 0),colspan=3))
            
            data=rk05[:,:]
            perc=np.array(df.bic_tp_area05)
            lb=list(nm)
            for i in range(len(nm)):
                lb[i]=str(nm[i]+":"+str(np.round(perc[i],1)) )
                
            
            
            m=Basemap(projection='cyl',
                      llcrnrlon=np.min(lon),urcrnrlon=np.max(lon), 
                        llcrnrlat=lat.min(),urcrnrlat=lat.max(), 
                        resolution='c')       
            cax=m.pcolormesh(lon,lat,data,shading='flat',cmap=cmapnm)
            plt.clim(cmin,cmax)
            m.drawparallels([0])
            m.drawmeridians([180])
            
            deg=np.arange(0,360,60)
            lblon=[None]*len(deg)
            for i,d in enumerate(deg): lblon[i]=lon2str(d)
            ax1.set_xticks(deg)
            ax1.set_xticklabels(lblon,fontsize=15)
            
            deg=np.arange(-60,61,30)
            lblat=[None]*len(deg)
            for i,d in enumerate(deg): lblat[i]=lat2str(d)
            ax1.set_yticks(deg)
            ax1.set_yticklabels(lblat,fontsize=15)
            
            cbar=plt.colorbar(cax,ticks=range(len(nm)), #shrink=0.7,
                              #label='Preferred Noise Model',fontsize=15,
                              orientation='vertical')
            cbar.ax.set_yticklabels(lb,fontsize=15) # horizontally oriented colorbar
            plt.ylabel('bic_tp\n',fontsize=18)
            ax1.text(0.08,0.90,' e ',weight='bold',fontsize=18,
                    horizontalalignment='left',
                    transform=ax1.transAxes)
            #plt.ylabel('Latitude\n',fontsize=15)
            #plt.xlabel('\nLongitude',fontsize=15)
            
            ##################
            
            
            #plt.show()
            #% %
            fig.savefig(pwd+'plot/{}.png'.format(name),
                    format='png',dpi=300,bbox_inches='tight')
            
            plt.close()
            
            #% %
            # find the best trend and unc
            trend_aic05=np.zeros((len(lat),len(lon)))
            trend_aic05.fill('nan')
            trend_bic05=np.zeros((len(lat),len(lon)))
            trend_bic05.fill('nan')
            trend_bic_tp05=np.zeros((len(lat),len(lon)))
            trend_bic_tp05.fill('nan')
            
            unc_aic05=np.zeros((len(lat),len(lon)))
            unc_aic05.fill('nan')
            unc_bic05=np.zeros((len(lat),len(lon)))
            unc_bic05.fill('nan')
            unc_bic_tp05=np.zeros((len(lat),len(lon)))
            unc_bic_tp05.fill('nan')
            for ilat,l in enumerate(lat):
                for ilon,l in enumerate(lon):
                        ind=np.array(rka05[ilat,ilon])
                        if np.isnan(ind):
                            trend_aic05[ilat,ilon]='nan'
                            unc_aic05[ilat,ilon]='nan'
                        else:
                            trend_aic05[ilat,ilon]=trend05[int(ind),ilat,ilon]
                            unc_aic05[ilat,ilon]=unc05[int(ind),ilat,ilon]
                        
                        ind=np.array(rkb05[ilat,ilon])
                        if np.isnan(ind):
                            trend_bic05[ilat,ilon]='nan'
                            unc_bic05[ilat,ilon]='nan'
                        else:
                            trend_bic05[ilat,ilon]=trend05[int(ind),ilat,ilon]
                            unc_bic05[ilat,ilon]=unc05[int(ind),ilat,ilon]
                        
                        ind=np.array(rk05[ilat,ilon])
                        if np.isnan(ind):
                            trend_bic_tp05[ilat,ilon]='nan'
                            unc_bic_tp05[ilat,ilon]='nan'
                        else:
                            trend_bic_tp05[ilat,ilon]=trend05[int(ind),ilat,ilon]
                            unc_bic_tp05[ilat,ilon]=unc05[int(ind),ilat,ilon]
                            
            best_trend_aic[:,:]=trend_aic05
            best_trend_bic[:,:]=trend_bic05
            best_trend_bic_tp[:,:]=trend_bic_tp05
            best_unc_aic[:,:]=unc_aic05
            best_unc_bic[:,:]=unc_bic05
            best_unc_bic_tp[:,:]=unc_bic_tp05        
            
            # add to dataframe
            ranks = np.zeros((3,len(lat),len(lon)))
            best_trend = np.zeros((3,len(lat),len(lon)))
            best_unc = np.zeros((3,len(lat),len(lon)))
            
            ranks[0,:,:]=rka
            ranks[1,:,:]=rkb
            ranks[2,:,:]=rk
            
            best_trend[0,:,:]=best_trend_aic
            best_trend[1,:,:]=best_trend_bic
            best_trend[2,:,:]=best_trend_bic_tp
            
            best_unc[0,:,:]=best_unc_aic
            best_unc[1,:,:]=best_unc_bic
            best_unc[2,:,:]=best_unc_bic_tp
            
            ds['ranks']=(('ic','lat','lon'),ranks)
            ds['best_trend']=(('ic','lat','lon'),best_trend)
            ds['best_unc']=(('ic','lat','lon'),best_unc)
            # ds = ds.assign_coords(ic=("ic", ['aic','bic','bic_tp']))
            
            # ds.to_netcdf(pwd+'ranking/{}.nc'.format(name))
    
        else:
            
            print('mean data (reg,nm)')
            #% % Rank noise models for each dataset
            reg=np.array(ds.reg)
            trend=np.array(ds.trend)
            uncert=np.array(ds.unc)
            aic=np.array(ds.aic)
            bic=np.array(ds.bic)
            bic_tp=np.array(ds.bic_tp)#np.nanmean([aic,bic],axis=0)
            nm=np.array(ds.nm)
            
            #% %
            rkb=np.full_like(np.zeros((len(reg))),np.nan)
            rka=np.full_like(rkb,np.nan)
            rk=np.full_like(rkb,np.nan)
            
            best_trend_aic=np.full_like(rkb,np.nan)
            best_trend_bic=np.full_like(rkb,np.nan)
            best_trend_bic_tp=np.full_like(rkb,np.nan)
            best_unc_aic=np.full_like(rkb,np.nan)
            best_unc_bic=np.full_like(rkb,np.nan)
            best_unc_bic_tp=np.full_like(rkb,np.nan)
            
            
            for ireg in range(len(reg)):
                print(ireg)
            
                targeta=aic[ireg,:] # All datasets and noise-models for this lat and lon
                targetb=bic[ireg,:]       
                targetab=bic_tp[ireg,:]
                
                bic2=np.zeros((len(nm)))
                aic2=np.zeros((len(nm)))
                bic_tp2=np.zeros((len(nm)))
            
                # for iname in np.arange(0,len([name])):
                iname = 0
            ##                AIC
                tg=targeta[:] # all noise models for this dataset
                logaic=np.zeros((tg.shape))
                for inm in np.arange(0,len(nm)):
                    logaic[inm]=np.exp((np.nanmin(tg)-tg[inm])/2)
                
                logaic 
                indaic=np.where(logaic>0.5)
                indaic=np.array(indaic)
                for inm in np.arange(0,len(nm)):
                    for ind in indaic:
                        for ind2 in ind:
                            if inm==ind2:        
                                aic2[inm]=1
            #                
                
                #BIC
                tg=targetb[:]
                tg
                logbic=np.zeros((tg.shape))
                for inm in np.arange(0,len(nm)):
                    logbic[inm]=np.exp((np.nanmin(tg)-tg[inm])/2)
                    # logbic[inm] = tg[inm] - np.nanmin(tg)
                
                logbic 
                indbic=np.where(logbic>0.5)
                # indbic=np.where(logbic<2)
                indbic=np.array(indbic)
                indbic
                for inm in np.arange(0,len(nm)):
                    for ind in indbic:
                        for ind2 in ind:
                            if inm==ind2:        
                                bic2[inm]=1
            #                AIC|BIC
                tg=targetab[:]
                tg
                logab=np.zeros((tg.shape))
                for inm in np.arange(0,len(nm)):
                    logab[inm]=np.exp((np.nanmin(tg)-tg[inm])/2)
                
                logab 
                indab=np.where(logab>0.5)
                indab=np.array(indab)
                indab
                for inm in np.arange(0,len(nm)):
                    for ind in indab:
                        for ind2 in ind:
                            if inm==ind2:        
                                bic_tp2[inm]=1
                                
            
                # indaic=np.sum(aic2,axis=0)                   
                # tmp=np.array(np.where(aic2==np.nanmax(indaic)))
                inda= np.where(aic2==np.max(aic2))
                tmp=np.array(np.where(targeta==targeta[inda].min()))
                if np.all(aic2==0):
                    rka[ireg]=0
                else:
                    rka[ireg]=tmp[0,0]
            
                # indbic=np.sum(bic2,axis=0)
                # tmp=np.array(np.where(bic2==np.nanmax(indbic)))
                indb= np.where(bic2==np.max(bic2))
                tmp=np.array(np.where(targetb==targetb[indb].min()))
                if np.all(bic2==0):
                    rkb[ireg]=0
                else:              
                    rkb[ireg]=tmp[0,0]
            
                # indab=np.sum(bic_tp2,axis=0)
                indab= np.where(bic_tp2==np.max(bic_tp2))
                tmp=np.array(np.where(targetab==targetab[indab].min()))
                if np.all(bic_tp2==0):
                    rk[ireg]=0
                else:
                    rk[ireg]=tmp[0,0]
                        
    
            #% %
            # rka05=np.array(rka)
            # rkb05=np.array(rkb)
            # rk05=np.array(rk)
            # trend05=np.array(trend)
            # unc05=np.array(uncert)
            nm_dic = {i:n for i,n in enumerate(nm)}
            df=pd.DataFrame(rkb,columns=['bic_rank'])
            df['bic_nm'] = [nm_dic[int(df['bic_rank'][i])] for i in range(len(rkb))]
            
            df['aic_rank']=rka
            df['aic_nm'] = [nm_dic[int(df['aic_rank'][i])] for i in range(len(rka))]
            
            df['bic_tp_rank']=rk
            df['bic_tp_nm'] = [nm_dic[int(df['bic_tp_rank'][i])] for i in range(len(rkb))]
            
            
            
            
            ax = plt.subplot(111, frame_on=False) # no visible frame
            ax.xaxis.set_visible(False)  # hide the x axis
            ax.yaxis.set_visible(False)  # hide the y axis
            
            table(ax, df, loc='center')  # where df is your data frame
            # #plt.show()
            plt.savefig(pwd+'plot/{}.png'.format(name))
            # fig.savefig(pwd+'plot/{}.png'.format(name),
            #         format='png',dpi=300,bbox_inches='tight')
            
            dfr=df
            
            #% %
            # find the best trend and unc
            trend_aic05=np.zeros((len(reg)))
            trend_aic05.fill('nan')
            trend_bic05=np.zeros((len(reg)))
            trend_bic05.fill('nan')
            trend_bic_tp05=np.zeros((len(reg)))
            trend_bic_tp05.fill('nan')
            
            unc_aic05=np.zeros((len(reg)))
            unc_aic05.fill('nan')
            unc_bic05=np.zeros((len(reg)))
            unc_bic05.fill('nan')
            unc_bic_tp05=np.zeros((len(reg)))
            unc_bic_tp05.fill('nan')
            for ireg in range(len(reg)):
                idx = int(df['aic_rank'][ireg])
                trend_aic05[ireg]=trend[ireg,idx]
                unc_aic05[ireg]=uncert[ireg,idx]
            
                idx = int(df['bic_rank'][ireg])
                trend_bic05[ireg]=trend[ireg,idx]    
                unc_bic05[ireg]=uncert[ireg,idx]
            
                idx = int(df['bic_tp_rank'][ireg])
                trend_bic_tp05[ireg]=trend[ireg,idx]   
                unc_bic_tp05[ireg]=uncert[ireg,idx]
            
            best_trend_aic[:]=trend_aic05
            best_trend_bic[:]=trend_bic05
            best_trend_bic_tp[:]=trend_bic_tp05
            best_unc_aic[:]=unc_aic05
            best_unc_bic[:]=unc_bic05
            best_unc_bic_tp[:]=unc_bic_tp05        
            
            # add to dataframe
            ranks = np.zeros((3,len(reg)))
            best_trend = np.zeros((3,len(reg)))
            best_unc = np.zeros((3,len(reg)))
            
            ranks[0,:]=rka
            ranks[1,:]=rkb
            ranks[2,:]=rk
            
            best_trend[0,:]=best_trend_aic
            best_trend[1,:]=best_trend_bic
            best_trend[2,:]=best_trend_bic_tp
            
            best_unc[0,:]=best_unc_aic
            best_unc[1,:]=best_unc_bic
            best_unc[2,:]=best_unc_bic_tp
            
            ds['ranks']=(('ic','reg'),ranks)
            ds['best_trend']=(('ic','reg'),best_trend)
            ds['best_unc']=(('ic','reg'),best_unc)
            
        ds = ds.assign_coords(ic=("ic", ['aic','bic','bic_tp']))
    
        ds.to_netcdf(pwd+'ranking/{}.nc'.format(name))        
