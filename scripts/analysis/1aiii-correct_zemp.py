#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 21:11:50 2022

@author: ccamargo
"""
import numpy as np
import xarray as xr
import pandas as pd
import os
from scipy.interpolate import interp1d
import pickle
# def main():
#     da = zemp_reg()
#     path = '/Volumes/LaCie_NIOZ/data/barystatic/revisions/input_data/'
#     da.to_netcdf(path+'zemp_month_reg.nc')
#     return 

import sys
sys.path.append("/Users/ccamargo/Documents/py_scripts/")
import utils_SLE_v2 as sle
import utils_SL as sl

import utils_hec as hec
# import os
# import cmocean as cmo
# from numba import jit
# import datetime as dt


def get_grid_area(grid):
    # given a grid that has dimensions: len(lat) x len(lon), compute the 
    #area in each grid cell, in km2
    # input: grid: numpy array
    
    earth_radius = 6371 # km
    # earth_radius = 6378137/1000# a more precise earth radius
    earth_diam = 2* earth_radius # diameter in km
    earth_circ = np.pi* earth_diam # earth's circunference in meters
    
    #% %
    dimlat,dimlon=grid.shape
    deltalat=180/dimlat
    deltalon=360/dimlon
    if deltalat==0.5:
        lat=np.arange(-89.875,90,deltalat)
        lon=np.arange(0.125,360,deltalon)
    else:
        lat=np.arange(-89.5,90,deltalat)
        lon=np.arange(0.5,360,deltalon)       

    #Transform from degrees to km:
    deltay=(earth_circ*deltalat)/360 #lat to km
    deltax=(earth_circ*np.cos(np.radians(lat))*deltalon)/360 #lon to km
    
    grid_area=np.array([deltax*deltay]*len(lon)).transpose()
    
    return grid_area

def glb_to_reg(value=1,mask = np.ones((180,360))):
    
    df=pd.DataFrame(mask.flatten(),columns=['mask'])
    df['area'] = get_grid_area(mask).flatten()
    
    df['regional_value'] = (np.full_like(mask,value).flatten() * df['mask'])/ len(mask[np.isfinite(mask)])# df['area']

    return np.array(df['regional_value']).reshape(mask.shape)


def load_dict(name, path ):
    with open(path + name + '.pkl', 'rb') as f:
        return pickle.load(f)
    
#% % open mask
def mask_gla():
    name = 'mask'
    path = '/Volumes/LaCie_NIOZ/data/barystatic/revisions/'
    
    dic = load_dict(name,path)
    # print(masks.keys())
    mask = dic['gla']['mask_regions']
    lat = dic['lat']
    lon=dic['lon']
    return mask,lat,lon

def year2mon(years,y):
    months = np.arange(
        np.datetime64('{}-01-01'.format(years[0])), 
        np.datetime64('{}-01-01'.format(years[-1]+1)),
        np.timedelta64(1,'M'),dtype='datetime64[M]')
    x = np.array(years)
    f = interp1d(x,y)
    xnew = np.linspace(years[0],years[-1],len(months))
    ynew = f(xnew)
    return xnew,ynew

#% %
def open_zemp():
    # dataset='GLA_ZMP'
    ymin=1962
    path = '/Volumes/LaCie_NIOZ/data/barystatic/original/GLA_Zemp2019-v1.0/'
    ext='.csv'
    flist=sorted([file for file in os.listdir(path) if file.endswith(ext) 
           and not file.startswith('.')])
    # sl.get_filelist('',)
    
    # global
    f = 'Zemp_etal_results_global.csv'
    df=pd.read_csv(path+f, header=18)
    df = df.set_index(df.Year)
    df = df.drop([year for year in df.index if year<ymin],axis =0 )
    years = np.array(df['Year'])
    y = np.array(df[' INT_SLE'])
    # yerr = np.array(df[' sig_Total_SLE'])
    # plt.plot(y);plt.plot(y+yerr);plt.plot(y-yerr)
    x = np.array(years)
    xnew,ynew = year2mon(x,y)
    # plt.plot(x, y, 'o', xnew, ynew, '-')
    flist.remove(f)
    
    # remove GRE and ANT:
    # flist.remove(f for f in ['Zemp_etal_results_region_19_ANT.csv',
    #                          'Zemp_etal_results_region_5_GRL.csv'])
    
    regs = [int(file.split('_')[4]) for file in flist]

    mchange_per_glacier_year=np.zeros((len(regs),len(years)))
    mchange_per_glacier_month=np.zeros((len(regs),len(xnew)))
    unc_per_glacier_year=np.zeros((len(regs),len(years)))
    unc_per_glacier_month=np.zeros((len(regs),len(xnew)))

    i = 0
    # plt.figure()
    for ireg, f in zip(regs,flist):
        df=pd.read_csv(path+f, header=26)
        df.head()
        df = df.set_index(df.Year)
        df = df.drop([year for year in df.index if year<ymin],axis =0 )
        x = np.array(df['Year'])
        # y = np.array(df[" AW_mwe"]) * 10**6# mwe = 10**3 kg/m2 = 10**6mm EWH
        mode='INT'
        unit='Gt'
        y = np.cumsum(np.array(df[" {}_{}".format(mode,unit)])/362.5 *10**5)# Gt to mm EWH
        # y = np.cumsum(np.array(df[" {}_{}".format(mode,unit)])/362.5)# Gt to mm SLE
        y = y-y[0]
        
        s = np.array(df[" sig_Total_{}".format(unit)])
        mchange_per_glacier_year[i] = np.array(y)
        _, ynew = year2mon(x,y)
        mchange_per_glacier_month[i] = ynew
        # plt.plot(x, y, 'o', xnew, ynew, '-')
        unc_per_glacier_year[i] = np.array(s)
        _, snew = year2mon(x,s)
        unc_per_glacier_month[i] = np.array(snew)
        # plt.plot(x, y+s, 'o', xnew, ynew+snew, '-')
        # plt.title(str(ireg))
        # plt.show()
        i=i+1
    return (xnew,regs,mchange_per_glacier_month,unc_per_glacier_month)

#%%
def zemp_reg():
    time,regs,mchange,unc = open_zemp()
    mask,lat,lon = mask_gla()
    ntime = len(time)
    nlat,nlon = mask.shape
    mchange2d = np.zeros((ntime,nlat*nlon))
    unc2d = np.full_like(mchange2d,0)
    
    for i,reg in enumerate(regs):
        if reg==19 or reg==5:
            # tmp = np.array(mask)
            # tmp[mask!=reg]=np.nan
            # tmp[np.isfinite(tmp)]=1
            pass
        else:
            tmp = np.array(mask)
            tmp[mask!=reg]=np.nan
            tmp[np.isfinite(tmp)]=1
            # plt.pcolor(tmp)
            for j in range(ntime):
                mchange2d[j,np.isfinite(tmp.flatten())] = glb_to_reg(mchange[i,j],tmp).flatten()[np.isfinite(tmp.flatten())]
                unc2d[j,np.isfinite(tmp.flatten())] = glb_to_reg(unc[i,j],tmp).flatten()[np.isfinite(tmp.flatten())]
    # Data is in mm of EWH
    
    # data = np.array(mchange2d).reshape(ntime,180,360)
    # tr, _ = sl.get_reg_trend(time,data,lat,lon)
    # # tr = np.array(data[-1])/55
    # # data = sle.height_to_EWH((tr).flatten()).reshape(180,360)
    # # data = np.array(tr*m*1000*10000)
    # _ = np.array(sle.run_SLE(tr,'zmp_tr',var='rsl')).reshape(180,360)
    #% % 
    da = xr.Dataset(data_vars={'mchange_2d_EWH':(('time','lat','lon'),mchange2d.reshape(ntime,nlat,nlon)),
                               'unc_2d_EWH':(('time','lat','lon'),unc2d.reshape(ntime,nlat,nlon)),
                               # 'unc_ewh':(('reg','time'),unc),
                               # 'mchange_SLE':(('reg','time'),mchange),
                               },
                    coords={'lat':lat,
                            'lon':lon,
                            'time':time,
                            'reg':regs})
    # data is in mm SLE
    for var in ['mchange_2d','unc_2d']:
        height = np.array(da[var+'_EWH'])
        data = np.full_like(height,0)
        for i in range(len(data)):
            data[i] = sle.EWH_to_height(height[i]).reshape(180,360)
        da[var+'_SLH'] = (('time','lat','lon'),data)
    da.attrs['SLH'] = 'mm of Sea-level Height'
    da.attrs['EWH'] = 'mm of Equivalent Water Height' 
    return da

#%%
def zemp_by_gla(da):
    mask,lat,lon = mask_gla()
    mask[mask==5] = np.nan
    mask[mask==19] = np.nan
    nglaciers = np.unique(mask[np.isfinite(mask)])
    ntime = len(da.time)
    table = np.zeros((len(nglaciers),ntime))
    table2 = np.zeros((len(nglaciers),ntime))
    u2 = np.zeros((len(nglaciers),ntime))
    u = np.zeros((len(nglaciers),ntime))
    
    # table2 = np.zeros((nglaciers,ntime))
    
    for i,g in enumerate(nglaciers):
        var = 'mchange_2d_EWH'
        m = np.array(mask)
        m[mask!=g] = np.nan
        m[np.isfinite(m)] = 1
        table[i] = np.nansum(np.array(da[var]*m),axis=(1,2))
        var = 'mchange_2d_SLH'
        # m = np.array(mask)
        # m[mask!=g] = np.nan
        # m[np.isfinite(m)] = 1
        table2[i] = np.nansum(np.array(da[var]*m),axis=(1,2))
        
        var = 'unc_2d_EWH'
        # m = np.array(mask)
        # m[mask!=g] = np.nan
        # m[np.isfinite(m)] = 1
        u[i] = np.nansum(np.array(da[var]*m),axis=(1,2))
        var = 'unc_2d_SLH'
        # m = np.array(mask)
        # m[mask!=g] = np.nan
        # m[np.isfinite(m)] = 1
        u2[i] = np.nansum(np.array(da[var]*m),axis=(1,2))
    
    ds = xr.Dataset(data_vars={'mchange_EWH':(('reg','time'),table),
                               'unc_EWH':(('reg','time'),u),
                               'mchange_SLH':(('reg','time'),table2),
                               'unc_SLH':(('reg','time'),u2),
                               'mask':(('lat','lon'),mask),
                               # 'unc_ewh':(('reg','time'),unc),
                               # 'mchange_SLE':(('reg','time'),mchange),
                               },
                    coords={'lat':lat,
                            'lon':lon,
                            'time':da.time,
                            'reg':nglaciers}) 
    ds.attrs['SLH'] = 'mm of Sea-level Height'
    ds.attrs['EWH'] = 'mm of Equivalent Water Height'
    
    return ds
#%% run
da = zemp_reg()
path = '/Volumes/LaCie_NIOZ/data/barystatic/revisions/input_data/'
da.to_netcdf(path+'zemp_month_reg.nc')
ds = zemp_by_gla(da)
ds.to_netcdf(path+'zemp_month_bygla.nc')

#%%
def run_hector():
    
    
    periods =[ (2005,2016),(1993,2018),
              (2003,2017)]
    for period in periods:
        t0,t1 = period
        print('{} to {}'.format(t0,t1-1))
    
        ifolder = 0    
        path_to_hec='/Volumes/LaCie_NIOZ/dump/'+str(ifolder)+'/'
        os.chdir(path_to_hec)
        
        #% % load dataset
        path = '/Volumes/LaCie_NIOZ/data/barystatic/revisions/input_data/'
        ds=xr.open_dataset(path+'zemp_month_bygla.nc')
        
    
        path_to_save ='/Volumes/LaCie_NIOZ/data/barystatic/revisions/hector/{}-{}/'.format(t0,t1-1)
        
        #% % means
        factor=10**16
        # Select time
        ds= ds.sel(time=slice(t0,t1))
        
        # test diff noise models:
        NM=['WN' ,
            'PL',
            'PLWN',
            'AR1',
            'AR5',
            'AR9',
            'ARF',
            'GGMWN',# 'FNWN','RWFNWN'
            ]
        
        
    
        time=np.array(ds.time)
        sp=30
            
        regions = np.array(ds['reg'])
        tr = np.zeros((len(regions),len(NM)))
        tr_err = np.full_like(tr, 0)
        aic = np.full_like(tr, 0)
        bic = np.full_like(tr, 0)
        bic_c = np.full_like(tr, 0)
        bic_tp = np.full_like(tr, 0)
        logL = np.full_like(tr, 0)
        N = np.full_like(tr, 0)
    
    
        dataset='GLA_ZMP'
        data = np.array(ds['mchange_SLH'])
        name = dataset
        os.chdir(path_to_hec)
        # if dataset=='GLA_ZMP_slc':
        for ireg,reg in enumerate(regions):
            print(reg)
            x=data[ireg,:] * factor
            # create a .mom file for it:
            hec.ts_to_mom(x[np.isfinite(x)],time[np.isfinite(x)],
                              sp=sp,path=str(path_to_hec+'raw_files/'),
                              name=str(name),ext='.mom')
            # loop over NM:
            for inm, n in enumerate(NM):
                print(n)
                if n=='ARF':
                    hec.create_estimatetrend_ctl_file(name,n,sp=sp,
                                                          GGM_1mphi = 6.9e-07,
                                                          seas=True,halfseas=True,
                                                          LikelihoodMethod='FullCov'
                                                          )
                else:
                    hec.create_estimatetrend_ctl_file(name,n,sp=sp,
                                                          GGM_1mphi = 6.9e-07,
                                                          seas=False,halfseas=False,
                                                          LikelihoodMethod='FullCov',
                                                          )
    
                # Run estimatetrend (hector)
                os.system('estimatetrend > estimatetrend.out')
    
                # get results:
                out=hec.get_results()
                #save results:
                if out[0]!=None:
                    tr[ireg,inm]=float(out[0])/factor
                    tr_err[ireg,inm]=float(out[1])/factor
                    N[ireg,inm]=out[2]
                    logL[ireg,inm]=out[3]
                    aic[ireg,inm]=out[4]
                    bic[ireg,inm]=out[5]
                    bic_tp[ireg,inm]=out[6]
                    bic_c[ireg,inm]=out[7]
                
                    
            
        
            da=xr.Dataset(data_vars={'trend':(('reg','nm'),tr),
                                     'unc':(('reg','nm'),tr_err),
                                     'aic':(('reg','nm'),aic),
                                     'bic':(('reg','nm'),bic),
                                     'bic_c':(('reg','nm'),bic_c),
                                     'bic_tp':(('reg','nm'),bic_tp),
                                     'logL':(('reg','nm'),logL),
                                     'N':(('reg','nm'),N),
                                    },
                                    coords={'nm':NM,
                                            'fname':dataset,
                                            'reg':regions })
        da.attrs['metadata']='barystatic sea-level trends in mm/y from {} to {}, obtained with Hector'.format(t0+1,t1)
        os.system(' cd {}'.format(path_to_save))
        da.to_netcdf(path_to_save+'{}_{}_NM'.format(dataset,len(NM))+'.nc')
        
run_hector()