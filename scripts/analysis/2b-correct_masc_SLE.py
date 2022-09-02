#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 14:55:45 2022

Separate mascons following reviewer's comments
Re-run SLE

@author: ccamargo
"""
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append("/Users/ccamargo/Documents/py_scripts/")
# import utils_SL as sl 
import utils_SLE_v2 as sle
#% % mask
import pickle
def load_dict(name, path):
    with open(path + name + '.pkl', 'rb') as f:
        return pickle.load(f)
    
def get_mask():
    name = 'mask'
    path = '/Volumes/LaCie_NIOZ/data/barystatic/revisions/'
    
    dic = load_dict(name,path)

    return dic

#% % epen trends

def open_trend(t0,t1):
    path = '/Volumes/LaCie_NIOZ/data/barystatic/revisions/results/{}-{}/'.format(t0,t1)
    file = 'source_trend_temporal_unc_{}-{}.nc'.format(t0,t1)

    ds = xr.open_dataset(path+file)
    return ds
def cor_zemp(ds):
    #% %
    dic = get_mask()
    mask = np.array(dic['gla']['mask_noIS'])
    mask[np.isfinite(mask)] = 1
    # mask = np.array(dic['mask'])
    # mask[mask==7] = 5
    # mask[mask!=5] = np.nan 
    # mask[np.isfinite(mask)] = 1
    names = np.array(ds.name)
    name = 'GLA_ZMP'
    idx = np.where(names == name)[0][0]
    # remove greenland:
    for var in ['best_trend','best_unc','ranks']:
        data = np.array(ds[var])
        for i in range(len(ds.ic)):
            data[idx,i,:,:] = np.array(data[idx,i,:,:]*mask)
        ds[var] = (ds[var].dims,data)
    for var in ['trend','unc']:
        data = np.array(ds[var])
        for i in range(len(ds.nm)):
            data[idx,i,:,:] = np.array(data[idx,i,:,:]*mask)
        ds[var] = (ds[var].dims,data)  
    
    return ds
def wouters():
    import pandas as pd
    df = pd.read_csv('/Users/ccamargo/Documents/PhD/Barystatic/Wouters2019 - Table1.csv', header=1)
    dic = get_mask()
    mask = dic['gla']['mask_regions']
    trend = np.full_like(mask,np.nan)
    unc = np.full_like(trend,np.nan)
    for igla in np.array(df['code']):
        area, tr, u = np.array(df.loc[df['code']==igla][['Area','Mass_change_trend','Mass_change_unc']])[0]
        # rel_area = np.array(area/len(mask[mask==igla]))
        rel_area = len(mask[mask==igla])
        trend[mask==igla] = (tr/362.5)/rel_area 
        unc[mask==igla] = (u/362.5)/rel_area 
    trend = sle.height_to_EWH(trend).reshape(180,360)
    unc = sle.height_to_EWH(unc).reshape(180,360)
        
    return trend, unc
    
def marzeion():
    import pandas as pd 
    df = pd.read_csv('/Users/ccamargo/Documents/PhD/Barystatic/Malles2021 - Table3.csv', header=1)
    dic = get_mask()
    mask = dic['gla']['mask_regions']
    trend = np.full_like(mask,np.nan)
    unc = np.full_like(trend,np.nan)
    for igla in np.array(df['code']):
        tr, u = np.array(df.loc[df['code']==igla][['trend','unc']])[0]
        # rel_area = np.array(area/len(mask[mask==igla]))
        rel_area = len(mask[mask==igla])
        trend[mask==igla] = (-tr)/rel_area 
        unc[mask==igla] = (u)/rel_area 
    trend = sle.height_to_EWH(trend).reshape(180,360)
    unc = sle.height_to_EWH(unc).reshape(180,360)
        
    return trend, unc
    
def get_external_glaciers():
    import pandas as pd
    path = '/Users/ccamargo/Documents/PhD/Barystatic/hugonnet2021/extended_data_tables/'
    file = 'ED_table1_formatted_lastcolumns.xlsx'
    df = pd.read_excel(path+file,header=1,
                       names = ['Region_ID','Region_name',
                                  'Area','Mean elevation change rate',
                                  '±', 'Mean elevation change uncertainty',
                                  'Mass change rate', '+/-', 'Mass change uncertainty'])
    df = df.drop([0]) # remove empty row
    df = df.drop([20,21]) # remove global values
    df = df.drop([5,19]) # remove AIS and GIS
    # df = pd.read_csv('/Users/ccamargo/Desktop/Wouters2019 - Table1.csv', header=1)
    dic = get_mask()
    mask = dic['gla']['mask_regions']
    trend = np.full_like(mask,np.nan)
    unc = np.full_like(trend,np.nan)
    for igla in np.array(df['Region_ID']):
        area, tr, u = np.array(df.loc[df['Region_ID']==igla][['Area','Mass change rate','Mass change uncertainty']])[0]
        # rel_area = np.array(area/len(mask[mask==igla]))
        rel_area = len(mask[mask==igla])
        trend[mask==igla] = (tr/362.5)/rel_area
        unc[mask==igla] = (u/362.5)/rel_area
    trend = sle.height_to_EWH(trend).reshape(180,360)
    unc = sle.height_to_EWH(unc).reshape(180,360)
        
    return trend, unc
def _cor_msc(ds):
    #% %
    dic = get_mask()
    mask = np.array(dic['mask'])
    m = np.array(dic['gla']['mask_regions'])
    m[m==5] = np.nan
    m[m==19] = np.nan
    m[np.isfinite(m)]=1
    m_full = np.array(m)
    m_full [mask!=2] = np.nan
    m_split = np.array(m)
    m_split[np.isfinite(m_full)] = np.nan
    
    
    m_split = np.array(mask)
    m_split[mask!=3]=np.nan 
    m_split[mask==3] = 1
    m_full = np.array(mask)
    m_full[mask!=2] = np.nan
    m_full[mask==2] = 1
    names = np.array(ds.name)
    name = 'GLA_ZMP'
    idx_zmp = np.where(names == name)[0][0]
    tr_cor,un_cor = get_external_glaciers()
    # tr_cor,un_cor = wouters()
    # name_msc = [name for name in names if name.split('_')[-1]=='CSR' or name.split('_')[-1]=='JPL']
    # name_msc = [name for name in name_msc if name.split('_'[]0)=='']
    # idxs = [ np.where(names == name)[0][0] for name in name_msc]
    for msc in ['JPL','CSR']:
        name = 'GLA_{}'.format(msc)
        idx_gla = np.where(names == name)[0][0]
        name = 'LWS_{}'.format(msc)
        idx_lws = np.where(names == name)[0][0]
        #% %
        for var in ['best_trend','best_unc','ranks']:
            if var.split('_')[-1]=='unc':
                data_cor = np.array(un_cor)
            else:
                data_cor = np.array(tr_cor)
            data = np.array(ds[var])
            for i in range(len(ds.ic)):
                # get external glacier estimate for mascons plit locations:
                glacier_ext = np.array(data_cor[idx_zmp][i]*m_split)

                # mascon total signal :
                mascon_gla = np.array(data[idx_gla][i])
                # glacier full = mascon total signal in these locations:
                glacier_full = np.array(mascon_gla * m_full)
                # final glacier is glacier full + glacier_ext
                glacier_cor = np.array(glacier_full)
                glacier_cor[np.isfinite(m_split)] = glacier_ext[np.isfinite(m_split)]
                glacier_cor[glacier_cor==0] = np.nan
                
                lws_split = np.full_like(mascon_gla,0)
                lws_split[np.isfinite(m_split)] = np.array(mascon_gla[np.isfinite(m_split)] - glacier_ext[np.isfinite(m_split)])
                lws_full = np.array(data[idx_lws][i])
                lws_full[np.isnan(lws_full)] = 0
                lws_cor = np.array(lws_full - lws_split)
                lws_cor[lws_cor==0] = np.nan
                
                
                data[idx_gla,i,:,:] = np.array(glacier_cor)
                data[idx_lws,i,:,:] = np.array(lws_cor)
                
                # old_gla = np.array(data[idx_gla][i])
                # new_gla = np.array(old_gla * m_full)
                # zmp_gla = np.array(data[idx_zmp][i]*m_split)
                # # gla_full[np.isnan(gla_full)] = 0
                # # new_gla = np.aray(gla_full)
                # new_gla[np.isfinite(m_split)] = np.array(zmp_gla[np.isfinite(m_split)])
                
                # old_lws = np.array(data[idx_lws][i])
                # old_gla = np.array(old_gla - zmp_gla * m_split)
                # old_gla[np.isnan(old_gla)] = 0
                # old_lws[np.isnan(old_lws)]= 0 
                # new_lws = np.array(old_lws - old_gla)
                # new_lws[new_lws==0] = np.nan
                
                # data[idx_gla,i,:,:] = np.array(new_gla)
                # data[idx_lws,i,:,:] = np.array(new_lws)
            ds[var] = (ds[var].dims,data)
        for var in ['trend','unc']:  
            data = np.array(ds[var])
            for i in range(len(ds.nm)):
                # old_gla = np.array(data[idx_gla][i])
                # new_gla= np.array(old_gla * m_full)
                # zmp_gla = np.array(data[idx_zmp][i]*m_split)
                # new_gla[np.isfinite(m_split)] = np.array(zmp_gla[np.isfinite(m_split)])

                # old_lws = np.array(data[idx_lws][i])
                # old_gla = np.array(old_gla - zmp_gla * m_split)
                # old_gla[np.isnan(old_gla)] = 0
                # old_lws[np.isnan(old_lws)]= 0 
                # new_lws = np.array(old_lws - old_gla)
                # new_lws[new_lws==0] = np.nan
                
                # data[idx_gla,i,:,:] = np.array(new_gla)
                # data[idx_lws,i,:,:] = np.array(new_lws)
                
                # get external glacier estimate for mascons plit locations:
                glacier_ext = np.array(data_cor[idx_zmp][i]*m_split)

                # mascon total signal :
                mascon_gla = np.array(data[idx_gla][i])
                # glacier full = mascon total signal in these locations:
                glacier_full = np.array(mascon_gla * m_full)
                # final glacier is glacier full + glacier_ext
                glacier_cor = np.array(glacier_full)
                glacier_cor[np.isfinite(m_split)] = glacier_ext[np.isfinite(m_split)]
                glacier_cor[glacier_cor==0] = np.nan
                
                lws_split = np.full_like(mascon_gla,0)
                lws_split[np.isfinite(m_split)] = np.array(mascon_gla[np.isfinite(m_split)] - glacier_ext[np.isfinite(m_split)])
                lws_full = np.array(data[idx_lws][i])
                lws_cor = np.array(lws_full - lws_split)
                lws_cor[lws_cor==0] = np.nan
                
                
                data[idx_gla,i,:,:] = np.array(glacier_cor)
                data[idx_lws,i,:,:] = np.array(lws_cor)
            ds[var] = (ds[var].dims,data)
    return ds

def cor_msc(ds):
    #% %
    dic = get_mask()
    mask = np.array(dic['mask'])
    m_split = np.array(mask)
    m_split[mask!=3]=np.nan
    m_split[mask==3] = 1
    m_full = np.array(mask)
    m_full[mask!=2] = np.nan
    m_full[mask==2] = 1
    names = np.array(ds.name)
    tr_w19,un_w19 = get_external_glaciers()
    # name = 'GLA_ZMP'
    # idx_zmp = np.where(names == name)[0][0]
    
    # name_msc = [name for name in names if name.split('_')[-1]=='CSR' or name.split('_')[-1]=='JPL']
    # name_msc = [name for name in name_msc if name.split('_'[]0)=='']
    # idxs = [ np.where(names == name)[0][0] for name in name_msc]
    for msc in ['JPL','CSR']:
        name = 'GLA_{}'.format(msc)
        idx_gla = np.where(names == name)[0][0]
        name = 'LWS_{}'.format(msc)
        idx_lws = np.where(names == name)[0][0]
        #% %
        for var in ['best_trend','best_unc',
                    # 'ranks'
                    ]:
            data = np.array(ds[var])
            if var.split('_')[-1]=='unc':
                data_w19 = np.array(un_w19)
            else:
                data_w19 = np.array(tr_w19)
            for i in range(len(ds.ic)):
                old_gla = np.array(data[idx_gla][i])
                new_gla = np.array(old_gla * m_full)
                # zmp_gla = np.array(data[idx_zmp][i]*m_split)
                gla_cor = np.array(data_w19*m_split)
                # gla_full[np.isnan(gla_full)] = 0
                # new_gla = np.aray(gla_full)
                new_gla[np.isfinite(m_split)] = np.array(gla_cor[np.isfinite(m_split)])
                
                old_lws = np.array(data[idx_lws][i])
                old_gla = np.array(old_gla - gla_cor * m_split)
                old_gla[np.isnan(old_gla)] = 0
                old_lws[np.isnan(old_lws)]= 0
                new_lws = np.array(old_lws + old_gla)
                new_lws[new_lws==0] = np.nan
                
                data[idx_gla,i,:,:] = np.array(new_gla)
                data[idx_lws,i,:,:] = np.array(new_lws)
            ds[var] = (ds[var].dims,data)
            
        for var in ['trend','unc']:
            data = np.array(ds[var])
            if var=='unc':
                data_w19 = np.array(un_w19)
            else:
                data_w19 = np.array(tr_w19)
                
            for i in range(len(ds.nm)):
                old_gla = np.array(data[idx_gla][i])
                new_gla= np.array(old_gla * m_full)
                # zmp_gla = np.array(data[idx_zmp][i]*m_split)
                gla_cor = np.array(data_w19 * m_split)
                new_gla[np.isfinite(m_split)] = np.array(gla_cor[np.isfinite(m_split)])

                old_lws = np.array(data[idx_lws][i])
                old_gla = np.array(old_gla - gla_cor * m_split)
                old_gla[np.isnan(old_gla)] = 0
                old_lws[np.isnan(old_lws)]= 0
                new_lws = np.array(old_lws + old_gla)
                new_lws[new_lws==0] = np.nan
                
                data[idx_gla,i,:,:] = np.array(new_gla)
                data[idx_lws,i,:,:] = np.array(new_lws)
            ds[var] = (ds[var].dims,data)
    return ds
    
def run_SLE(ds,path_save,
            name_save,var='asl'):
    idx=-1
    # var = 'asl'
    slf_tr = np.full_like(np.zeros((len(ds.name),len(ds.lat),len(ds.lon))), 0 )
    slf_unc = np.full_like(slf_tr,0)
    # X=np.array(ds.lon); Y=np.array(-ds.lat)
    for iname,name in enumerate(np.array(ds.name)):
        slc = np.array(ds['best_trend'][iname][idx])
        unc = np.array(ds['best_unc'][iname][idx])
        
        if name.split('_')[-1]=='IMB':
            slc=-slc
        
        slf_tr[iname] = np.array(sle.run_SLE(slc,name+'_tr',var=var)).reshape(180,360)
        
        slf_unc[iname] = np.abs(np.array(sle.run_SLE(unc,name+'_unc',var=var)).reshape(180,360))
    
    da = xr.Dataset(data_vars = {'SLF_trend':(('name','lat','lon'),slf_tr),
                                 'SLF_unc':(('name','lat','lon'),slf_unc),
                                 },
                    coords = {'name':ds.name,
                              'lat':-ds.lat,
                              'lon':ds.lon}
                    )
    da.to_netcdf(path_save+name_save+'.nc')


def plot(ds):
    from cartopy import crs as ccrs , feature as cfeature
    import matplotlib.pyplot as plt
    import string
    # user defined functions
    import sys
    sys.path.append("/Users/ccamargo/Documents/py_scripts/OM/")
    import utils_OM as om
    import cmocean as cm
    from cmcrameri import cm as cmf
    from matplotlib.gridspec import GridSpec

    dpi=300
    wi=20;hi=15
    dimlat=180;dimlon=360
    fontsize=25
    ticksize=20
    
    letters = list(string.ascii_lowercase)
    landcolor='darkgrey'
    nrows=4
    ncols=4
    fig = plt.figure(figsize=(15,10), facecolor='w')
    interval=0.1
    
    gs = GridSpec(nrows, ncols)
    irow = 0
    
    # alpha=0.9
    lat=np.array(ds.lat)
    lon=np.array(ds.lon)
    # regions = ['AIS','GIS','GLA','LWS']
    # combos= ['JPL','UCI+WGP']
    panels=['Trend', 'Uncertainty']
    name_list = ['AIS_JPL','GIS_JPL',
                 'AIS_UCI','GIS_UCI',
                 'GLA_JPL','LWS_JPL',
                 'GLA_WGP','LWS_WGP']
    # ncols=len(panels)
    iletter=0
    for iname, name in enumerate(name_list):
        da=ds.sel(name=name)
        
        proj = ccrs.Robinson()
        for ipanel, panel in enumerate(panels):
            # print(icol)
            # print(col)
            ax1 = plt.subplot(gs[ipanel+irow], 
                                    projection=proj)
            ax1.coastlines(resolution='110m', zorder=3,color=landcolor)
            ax1.add_feature(cfeature.LAND,color=landcolor,# alpha=0.5,
                        zorder=3)
            ax1.set_global() 
    
            if irow==0 or irow==2:
                plt.title(#'{}\n{}'.format(panel,name),
                          '{}\n({}). {}'.format(panel,letters[iletter],name),
                          fontsize=20)
            else:
                plt.title('\n({}). '.format(letters[iletter])+name,fontsize=20)
            iletter=iletter+1
                
            ic_idx=-1
            if ipanel == 0 or ipanel==2:
                Y=np.array(lat);X=np.array(lon);X[0]=0;X[-1]=360
                clim=1
                cmap=cmf.roma_r
                data = np.array(da['SLF_trend'][:,:])  
    
    
                mu,glb,mask=om.reg_to_glb(data,Y,X)
                glb=np.round(glb,3)
                cmin=-clim;cmax=clim
                cs=ax1.contour(X,Y,data,levels=[glb],
                        # vmin=-0.6,
                        # vmax=0.6,
                    transform = ccrs.PlateCarree(),
                    #cmap='coolwarm',#extend='both'
                    colors=('black',),linestyles=('--',),linewidths=(2,)
                    )
                lv=np.arange(cmin,cmax+interval,interval)
                csf=plt.contourf(X,Y,data,levels=lv,
                          transform = ccrs.PlateCarree(),cmap=cmap)
                ax1.clabel(cs,cs.levels,fmt='%5.2f',colors='k',fontsize=12)
    
                plt.pcolormesh(X,Y,data,
                        vmin=cmin,vmax=cmax,
                        zorder=0,
                        transform = ccrs.PlateCarree(),cmap=cmap)
    
            else:
                data = np.abs(da['SLF_unc'][:,:])
    
                clim=1
                cmap=cmf.roma_r
                mu,glb,mask=om.reg_to_glb(data,Y,X)
                glb=np.round(glb,3)
                cmin=-clim;cmax=clim
               
                cs=ax1.contour(X,Y,data,levels=[glb],
                        # vmin=-0.6,
                        # vmax=0.6,
                    transform = ccrs.PlateCarree(),
                    #cmap='coolwarm',#extend='both'
                    colors=('black',),linestyles=('--',),linewidths=(2,)
                    )
                lv=np.arange(cmin,cmax+interval,interval)
                csf=plt.contourf(X,Y,data,levels=lv,
                          transform = ccrs.PlateCarree(),cmap=cmap)
                ax1.clabel(cs,cs.levels,fmt='%5.2f',colors='k',fontsize=12)
    
                cp2=plt.pcolormesh(X,Y,data,
                        vmin=cmin,vmax=cmax,
                        zorder=0,
                        transform = ccrs.PlateCarree(),cmap=cmap)
        irow = irow+len(panels)
               

    # xmin, xmax = ax1.get_xprop()
    # ymin, ymax = ax1.get_yprop()
    # y2x_ratio = (ymax-ymin)/(xmax-xmin) * nrows/ncols
    # fig.set_figheight(wi * y2x_ratio + 5)
    
    # fig.subplots_adjust(right=0.8)
    
    plt.tight_layout()
    # # fig.subplots_adjust(right=0.8)
    cbar_ax2 = fig.add_axes([0.025, -0.05, 0.95, 0.04])
    cbar2=plt.colorbar(csf, cax=cbar_ax2,orientation='horizontal')
    cbar2.set_label(label='mm/yr',size=ticksize, family='serif')    
    cbar2.ax.tick_params(labelsize=ticksize) 
    
    plt.show()

def plot2(ds):
    from cartopy import crs as ccrs , feature as cfeature
    # import matplotlib.pyplot as plt
    import string
    # user defined functions
    import sys
    sys.path.append("/Users/ccamargo/Documents/py_scripts/OM/")
    import utils_OM as om
    # import cmocean as cm
    from cmcrameri import cm as cmf
    from matplotlib.gridspec import GridSpec


    ticksize=20
    
    letters = list(string.ascii_lowercase)
    landcolor='darkgrey'
    nrows=4
    ncols=4
    fig = plt.figure(figsize=(15,10), facecolor='w')
    interval=0.1
    
    gs = GridSpec(nrows, ncols)
    irow = 0
    
    # alpha=0.9
    lat=np.array(ds.lat)
    lon=np.array(ds.lon)
    # regions = ['AIS','GIS','GLA','LWS']
    # combos= ['JPL','UCI+WGP']
    panels=['Trend', 'Uncertainty']
    name_list = ['AIS_CSR','GIS_CSR',
             'AIS_IMB','GIS_IMB',
             'GLA_CSR','LWS_CSR',
             'GLA_ZMP','LWS_GWB']
    # ncols=len(panels)
    iletter=0
    for iname, name in enumerate(name_list):
        da=ds.sel(name=name)
        
        proj = ccrs.Robinson()
        for ipanel, panel in enumerate(panels):
            # print(icol)
            # print(col)
            ax1 = plt.subplot(gs[ipanel+irow], 
                                    projection=proj)
            ax1.coastlines(resolution='110m', zorder=3,color=landcolor)
            ax1.add_feature(cfeature.LAND,color=landcolor,# alpha=0.5,
                        zorder=3)
            ax1.set_global() 
    
            if irow==0 or irow==2:
                plt.title(#'{}\n{}'.format(panel,name),
                          '{}\n({}). {}'.format(panel,letters[iletter],name),
                          fontsize=20)
            else:
                plt.title('\n({}). '.format(letters[iletter])+name,fontsize=20)
            iletter=iletter+1
                
            ic_idx=-1
            if ipanel == 0 or ipanel==2:
                Y=np.array(lat);X=np.array(lon);X[0]=0;X[-1]=360
                clim=1
                cmap=cmf.roma_r
                data = np.array(da['SLF_trend'][:,:])  
    
    
                mu,glb,mask=om.reg_to_glb(data,Y,X)
                glb=np.round(glb,3)
                cmin=-clim;cmax=clim
                cs=ax1.contour(X,Y,data,levels=[glb],
                        # vmin=-0.6,
                        # vmax=0.6,
                    transform = ccrs.PlateCarree(),
                    #cmap='coolwarm',#extend='both'
                    colors=('black',),linestyles=('--',),linewidths=(2,)
                    )
                lv=np.arange(cmin,cmax+interval,interval)
                csf=plt.contourf(X,Y,data,levels=lv,
                          transform = ccrs.PlateCarree(),cmap=cmap)
                ax1.clabel(cs,cs.levels,fmt='%5.2f',colors='k',fontsize=12)
    
                plt.pcolormesh(X,Y,data,
                        vmin=cmin,vmax=cmax,
                        zorder=0,
                        transform = ccrs.PlateCarree(),cmap=cmap)
    
            else:
                data = np.abs(da['SLF_unc'][:,:])
    
                clim=1
                cmap=cmf.roma_r
                mu,glb,mask=om.reg_to_glb(data,Y,X)
                glb=np.round(glb,3)
                cmin=-clim;cmax=clim
               
                cs=ax1.contour(X,Y,data,levels=[glb],
                        # vmin=-0.6,
                        # vmax=0.6,
                    transform = ccrs.PlateCarree(),
                    #cmap='coolwarm',#extend='both'
                    colors=('black',),linestyles=('--',),linewidths=(2,)
                    )
                lv=np.arange(cmin,cmax+interval,interval)
                csf=plt.contourf(X,Y,data,levels=lv,
                          transform = ccrs.PlateCarree(),cmap=cmap)
                ax1.clabel(cs,cs.levels,fmt='%5.2f',colors='k',fontsize=12)
    
                cp2=plt.pcolormesh(X,Y,data,
                        vmin=cmin,vmax=cmax,
                        zorder=0,
                        transform = ccrs.PlateCarree(),cmap=cmap)
        irow = irow+len(panels)
               

    # xmin, xmax = ax1.get_xprop()
    # ymin, ymax = ax1.get_yprop()
    # y2x_ratio = (ymax-ymin)/(xmax-xmin) * nrows/ncols
    # fig.set_figheight(wi * y2x_ratio + 5)
    
    # fig.subplots_adjust(right=0.8)
    
    plt.tight_layout()
    # # fig.subplots_adjust(right=0.8)
    cbar_ax2 = fig.add_axes([0.025, -0.05, 0.95, 0.04])
    cbar2=plt.colorbar(csf, cax=cbar_ax2,orientation='horizontal')
    cbar2.set_label(label='mm/yr',size=ticksize, family='serif')    
    cbar2.ax.tick_params(labelsize=ticksize) 
    
    plt.show()
#%% run 
periods = [
            (2003,2016), 
             #(1993,2017)
            ]
for period in periods:
    t0,t1 = period
    # var='rsl'
    
    ds = open_trend(t0,t1)
    # ds = cor_zemp(ds)
    if t0>2002:
        ds = cor_msc(ds)
    path_save = '/Volumes/LaCie_NIOZ/data/barystatic/revisions/results_final/{}-{}/'.format(t0,t1)
    ds.to_netcdf(path_save+'source_trend_unc_{}-{}.nc'.format(t0,t1))
    for var in ['rsl',
                # 'asl'
                ]:
        name_save = 'SLF_{}_{}-{}'.format(var,t0,t1)
    
        run_SLE(ds,path_save,name_save,var=var)
        da = xr.open_dataset(path_save+name_save+'.nc')
    #plot(da)
    # plot2(da)
#%% test 
# t0,t1  = (2003,2016)
# ds = open_trend(t0,t1)
# ds = cor_msc(ds)
# var='rsl'
# idx=-1
# # var = 'asl'
# names = np.array(ds.name)
# for iname in [
#                8,9,
#                12,
#                14
#               ]:
#     name = names[iname]
#     slc = np.array(ds['best_trend'][iname][idx])
#     _ = np.array(sle.run_SLE(slc,name+'_tr',var=var)).reshape(180,360)
    
    
