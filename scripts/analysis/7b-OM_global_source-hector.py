#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 21 10:44:18 2021

Hector on OM global mean time series 
also OLS 

find best NM and make plot

@author: ccamargo
"""


#% % libraries
# hector_reg_bary_9
import numpy as np
import xarray as xr
import sys
sys.path.append("/Users/ccamargo/Documents/py_scripts/")
import utils_SL as sl 
import utils_hec as hec
import os
# import cmocean as cmo
# from numba import jit
import datetime as dt

import utils_SLE_v2 as sle

import matplotlib.pyplot as plt
#%%
def find_best_val(tr,unc,idx):
    """
    tr,unc and idx are 1D arrays with the length of the used noise models.
    tr are the trends for each NM, unc the unc for each NM, 
    and idx is the AIC, BIC or BICtp ranking for each noise model
    """
    target = idx # all NM for this time series
    len_nm= len(idx)
    logidx= np.array([np.exp((np.nanmin(target)-target[inm])/2) for inm in range(0,len_nm)])
    
    threshold = 0.5
    if np.any(logidx>threshold):
        ind=np.where(logidx>0.5)[0]
        if len(ind)>1:
            ind=np.where(tr==np.min(tr[ind]))[0][0]
        else:
            ind=ind[0]
    else:
        ind=np.where(np.nanmin(target))[0][0]
    
    best_tr=tr[ind]
    best_unc=unc[ind]
    
    return best_tr,best_unc,ind
#%% 
col_dict={1:"black", # WN
          2:"palegoldenrod", # PL
          3:"lightpink", # PLWN
          4:"orange", # AR1
          5:"teal", # Ar5
          6:"darkmagenta", # AR9
          7:"skyblue", # ARf
          8:"crimson", # GGM
          9:"grey" # OLS
          }
#%%
periods =[ 
            (2005,2016),
            (1993,2018),
            (1993,2017),
            (2003,2017)
          ]
period=periods[0]
for period in periods:
    t0=period[0]
    t1=period[1]
    print('{} to {}'.format(t0,t1-1))

    ifolder = 0
    # t0=2005
    # t1=2016 # until december of previous year
    # pwd='/export/lv1/user/ccamargo/OM/'
    pwd='/Volumes/LaCie_NIOZ/data/barystatic/global/'
    # dataset='means_EWH.nc'
    dataset='global_mean_timeseries_v2.nc'
    #dataset='3-yearly_SLF_monthly_1993-2020_180x360.nc'



    #% % set path to hector
    
    path_to_hec='/Volumes/LaCie_NIOZ/dump/'+str(ifolder)+'/'
    os.chdir(path_to_hec)
    #% % load dataset
    path=pwd
    ds=xr.open_dataset(path+dataset)
    
    # path_to_save ='//export/lv1/user/ccamargo/OM/hector/'+str(t0)+'-2016/'
    # path_to_save ='//export/lv1/user/ccamargo/OM/hector/{}-{}/'.format(t0,t1-1)
    # os.system(' cd {}'.format(path_to_save))
    
    # path_to_save ='/Volumes/LaCie_NIOZ/data/barystatic/hector/source/EWH/{}-{}/'.format(t0,t1-1)
    path_to_save =pwd+'hector/'
    
    #% % means
    factor=10**16
    # Select time
    ds= ds.sel(years=slice(t0,t1-1))
    ds['months'],_=sl.get_dec_time(np.array(ds.months))
    ds= ds.sel(months=slice(t0,t1))
    
    # get variables
    variables = list(ds.keys())
    datasets=[name for name in np.array(ds.name_monthly)]
    if t0<2002:
        datasets = [name for name in np.array(ds.name_monthly) 
                    if not name.endswith('CSR') and not name.endswith('JPL')]
    datasets.extend(name for name in np.array(ds.name_yearly))
    # test diff noise models:
    NM=['WN' ,
        'PL',
        'PLWN',
        'AR1',
        'AR5',
        'AR9',
        'ARF',
        'GGMWN',# 'FNWN','RWFNWN'
        'OLS'
        ]
    
    tr = np.zeros((len(datasets),len(NM)))
    tr_err = np.full_like(tr, 0)
    aic = np.full_like(tr, 999999)
    bic = np.full_like(tr, 999999)
    bic_c = np.full_like(tr, 999999)
    bic_tp = np.full_like(tr, 999999)
    logL = np.full_like(tr, 999999)
    N = np.full_like(tr, 999999)
    
    best_tr=np.zeros((len(datasets),3))
    best_unc=np.full_like(best_tr,0)
    best_nm=np.full_like(best_tr,0)
    
    # iname=13;name=datasets[iname]
    for iname, name in enumerate(datasets):
        print(name)
        # % %
        if name in np.array(ds.name_yearly):
            time=np.array(ds.years)
            sp=365
            idx = np.where(name==np.array(ds.name_yearly))[0][0]
            data = np.array(ds.OM_yearly_slc[:,idx])
            seas=False;halfseas=False
            
        elif name in np.array(ds.name_monthly):
            time=np.array(ds.months)
            sp=30    
            idx = np.where(name==np.array(ds.name_monthly))[0][0]
            data = np.array(ds.OM_monthly_slc[:,idx])
            seas=True;halfseas=True

        
        x=data * factor
        # create a .mom file for it:
        hec.ts_to_mom(x[np.isfinite(x)],time[np.isfinite(x)],
                          sp=sp,path=str(path_to_hec+'raw_files/'),
                          name=str(name),ext='.mom')
        # loop over NM:
        inm=-1;n=NM[inm]
        for inm, n in enumerate(NM):
            if n =='OLS':
                out = sl.get_OLS_trend(time[np.isfinite(x)], x[np.isfinite(x)])
            elif n=='ARF':
                hec.create_estimatetrend_ctl_file(name,n,sp=sp,
                                                      GGM_1mphi = 6.9e-07,
                                                      seas=True,halfseas=True,
                                                      LikelihoodMethod='FullCov'
                                                      )
                # Run estimatetrend (hector)
                os.system('estimatetrend > estimatetrend.out')
    
                # get results:
                out=hec.get_results()


            else:
                hec.create_estimatetrend_ctl_file(name,n,sp=sp,
                                                  GGM_1mphi = 6.9e-07,
                                                      seas=seas,halfseas=halfseas,
                                                      LikelihoodMethod='FullCov',
                                                      )
                # Run estimatetrend (hector)
                os.system('estimatetrend > estimatetrend.out')
    
                # get results:
                out=hec.get_results()

            #save results:
            
            if n=='OLS': # OLS
                tr[iname,inm]=out[0]/factor
                tr_err[iname,inm]=out[1]/factor
            elif out[0]!=None:
                    tr[iname,inm]=float(out[0])/factor
                    tr_err[iname,inm]=float(out[1])/factor
                    N[iname,inm]=out[2]
                    logL[iname,inm]=out[3]
                    aic[iname,inm]=out[4]
                    bic[iname,inm]=out[5]
                    bic_tp[iname,inm]=out[6]
                    bic_c[iname,inm]=out[7]
        # end loop over NM
                    
        #% % Quick QC
        std = np.nanstd(tr[iname,:]) + np.nanmax(tr_err[iname,:])
        mu = np.nanmean(tr[iname,:])
        if np.any(tr[iname,:]> mu +std) or np.any(tr[iname,:]< mu-std):
            idx=np.array(np.where((tr[iname,:] > mu+std) | (tr[iname,:] < mu-std)))[0]
            for ix in idx:
                # print('QC removed {}'.format(nm[ix]))
                tr[iname,ix]=np.nan;
                tr_err[iname,ix]=np.nan
                # aic[jname,ix,ilat,ilon]=np.nan;
                # bic[jname,ix,ilat,ilon]=np.nan;
                # logL[jname,ix,ilat,ilon]=np.nan
                aic[iname,ix]=np.nanmax(aic[iname,:])
                bic[iname,ix]=np.nanmax(bic[iname,:])
                bic_tp[iname,ix]=np.nanmax(bic_tp[iname,:])
                logL[iname,ix]=np.nanmax(logL[iname,:])
                        
        # find best trend:
            ## for AIC:
        best_tr[iname,0],best_unc[iname,0],best_nm[iname,0] = find_best_val(tr[iname,:], tr_err[iname,:], aic[iname,:])
            ## for BIC:
        best_tr[iname,1],best_unc[iname,1],best_nm[iname,1] = find_best_val(tr[iname,:], tr_err[iname,:], bic[iname,:])
            ## for BIC_tp:
        best_tr[iname,2],best_unc[iname,2],best_nm[iname,2] = find_best_val(tr[iname,:], tr_err[iname,:], bic_tp[iname,:])
        # make dataset:
    da=xr.Dataset(data_vars={'trend':(('name','nm'),tr),
                                     'unc':(('name','nm'),tr_err),
                                     'aic':(('name','nm'),aic),
                                     'bic':(('name','nm'),bic),
                                     'bic_c':(('name','nm'),bic_c),
                                     'bic_tp':(('name','nm'),bic_tp),
                                     'logL':(('name','nm'),logL),
                                     'N':(('name','nm'),N),
                                     'best_tr':(('name','ic_idx'),best_tr),
                                     'best_unc':(('name','ic_idx'),best_unc),
                                     'best_nm':(('name','ic_idx'),best_nm)
                                    },
                                    coords={'nm':NM,
                                            'name':datasets,
                                            'ic_idx':['aic','bic','bic_tp'],
                                            })
    #% %
    da=da.sortby(da.name)
    
    # % %  PLOT
    tr=np.array(da.trend)
    unc=np.array(da.unc)
    nm=np.array(da.nm)
    name=np.array(da.name)
    idx=-1
    best_tr=np.array(da.best_tr[:,idx])
    best_unc=np.array(da.best_unc[:,idx])
    best_nm=np.array(da.best_nm[:,idx])
    #% % plot
    ind=np.arange(0,len(nm))
    x=np.arange(0,len(name))
    syb=['o','v',"^",'s','p','*','D','X','X','X']
    # cor =
    fig=plt.figure(figsize=(10,10),dpi=100,facecolor='w')
    ax=plt.subplot(3,1,1)
    for i in ind:
        y=tr[:,i]
        plt.scatter(x,np.abs(y),color=col_dict[i+1],label=nm[i],marker=syb[i],s=80)
    # plt.plot(x,np.repeat(0,len(x)),color='red',linestyle='--')
    plt.ylabel('Trend (mm/yr)',fontsize=15)
    plt.title('Barystatic Sea-level\n{}-{}'.format(t0,t1),fontsize=20)
    # plt.ylim(0,1)
    # plt.xticks(x,'',rotation=75,fontsize=15)
    plt.yticks(fontsize=15)
    plt.xticks(x,name,rotation=75,fontsize=15)
    plt.grid()
    plt.tight_layout()
    
    ax=plt.subplot(3,1,2)
    for i in ind:
        y=unc[:,i]
        plt.scatter(x,np.abs(y),color=col_dict[i+1],label=nm[i],marker=syb[i],s=80)
    ax.legend(bbox_to_anchor=(1.2, 1.5))
    # plt.ylim(0,0.55)
    # plt.plot(x,np.repeat(0,len(x)),color='red',linestyle='--')
    plt.grid()
    # plt.legend(loc='upper left',fontsize=15,ncol=2)
    plt.xticks(x,name,rotation=75,fontsize=15)
    plt.ylabel('Uncertainty (mm/yr)',fontsize=15)
    plt.yticks(fontsize=15)
    plt.ylim(-0.1,2)
    #  

    ax=plt.subplot(3,1,3)
    # plt.scatter(x+0,best_tr,#marker=syb,
    #             s=100,
    #             # color=col_dict[best_nm+1]
    #             )
    for i in x:
        print('{}: {}±{}'.format(name[i],np.round(best_tr[i],3),np.round(best_unc[i],3)))
        y=best_tr[i]
        plt.scatter(i,np.abs(y),color=col_dict[best_nm[i]+1],
                    # label=nm[best_nm[i]],
                    marker=syb[int(best_nm[i])],
                    s=150
                    )
        yerr=best_unc[i]
        plt.errorbar(i,np.abs(y),yerr=yerr,
                 # color=df['cor'],
                 capsize=3,capthick=2,ecolor='black',lw=2,fmt='none')
    # plt.plot(x,np.repeat(0,len(x)),color='red',linestyle='--')
    plt.grid()
    # plt.ylim(-1,1)
    # plt.legend(loc='upper left',fontsize=15,ncol=2)
    plt.xticks(x,name,rotation=75,fontsize=15)
    plt.ylabel('Best trend ± unc (mm/yr)',fontsize=15)
    plt.yticks(fontsize=15)
    #% % 
    fig.savefig(path_to_save+'plot/'+'ALL-ocean_mean_{}_NM-{}-{}'.format(len(NM),t0,t1-1)+'.png',
                format='png')
    plt.close()
    # SAVE
    da.attrs['metadata']='ocean mean barystatic sea-level trends in mm/y from {} to {}, obtained with Hector'.format(t0,t1-1)
    # os.system(' cd {}'.format(path_to_save))
    da.to_netcdf(path_to_save+'ALL-ocean_mean_{}_NM-{}-{}'.format(len(NM),t0,t1-1)+'_v2.nc')

