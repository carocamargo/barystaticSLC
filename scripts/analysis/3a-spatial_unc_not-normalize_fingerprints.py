#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 16:27:50 2022

Test Spatial Uncertainty without normalizing 

@author: ccamargo
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  7 17:29:50 2022

Spatial Uncertainty

Normalized Fingerprints

@author: ccamargo
"""


import sys
sys.path.append("/Users/ccamargo/Documents/py_scripts/")
import utils_unc as unc
import utils_SL as sl
import utils_SLE_v2 as sle

import pandas as pd
import xarray as xr
import numpy as np 
# import cmocean as cm
# from cmcrameri import cm as cmf
# import matplotlib.pyplot as plt

# %% # open masks:
# ANT 
path='/Volumes/LaCie_NIOZ/data/barystatic/masks/rignot_basins/ANT_Basins_IMBIE2_v1.6/'
ds=xr.open_dataset(path+'final_mask.nc')
ds
lon=np.array(ds.lon)
lon=sl.from_180_to_360(lon)
ds = ds.assign_coords(lon=lon)
ds = ds.sortby('lon')
lon=np.array(ds.lon)
lat=np.array(ds.lat)
# plot_ant(ds.lon,ds.lat,ds.mask,cmin=0,cmax=19,title='Rignot Drainage Basins',cmap='tab20')
mask=np.array(ds.mask)
codes=np.unique(mask[np.isfinite(mask)])
maskant_basins=np.array(ds.mask)
#% %
maskant_regions = np.array(maskant_basins)
wais=(1,2,9,14,18)
for i in wais:
    maskant_regions[np.where(maskant_basins==i)]=1
eais=(3,4,5,6,7,8,10,11,12,13)
for i in eais:
    maskant_regions[np.where(maskant_basins==i)]=2
ap=(15,16,17)
for i in ap:
    maskant_regions[np.where(maskant_basins==i)]=3    

# plot_ant(ds.lon,ds.lat,maskant_regions,cmin=0,cmax=3,
#          title='Rignot Drainage Basins per Region',cmap='tab10')

maskant = np.array(mask)
maskant[np.isfinite(mask)]=1
# plot_ant(ds.lon,ds.lat,maskant,cmin=0,cmax=1,title='AIS mask')
ngrid = np.array(maskant)
ngrid=np.hstack(ngrid)
j=0
for i in range(len(ngrid)):
    if np.isfinite(ngrid[i]):
        ngrid[i]=j
        j=j+1
ngrid_ant=ngrid.reshape(maskant.shape)

#%% GRE mask
ds=xr.open_dataset('/Volumes/LaCie_NIOZ/data/barystatic/masks/barystatic_mask2.nc')
maskgre_basins=np.array(ds.gre)
lat=np.array(ds.lat);lon=np.array(ds.lon)
dimlat=len(lat);dimlon=len(lon)
# unc.plot_gre(lon,lat,maskgre_basins,cmin=0,cmax=8,
#              cmap='tab10',title='Rinot/Mouginot GIS mask')

# We want only the dynamic regions NW (0), CW(2) and SE (6)
maskgre_dyn=np.copy(maskgre_basins)
for i in [1,3,4,5]:
    # print(i)
    maskgre_dyn[np.where(maskgre_basins==i)]=10
maskgre_dyn[maskgre_dyn!=10]=0
maskgre_dyn[maskgre_dyn==10]=1
maskgre_dyn[np.isnan(maskgre_basins)]=np.nan
# unc.plot_gre(lon,lat,maskgre_dyn,cmin=0,cmax=1,cmap='tab10',title='Dynamic GIS mask')

maskgre=np.array(maskgre_basins)
maskgre[np.isfinite(maskgre)]=1
# unc.plot_gre(lon,lat,maskgre,title='GIS mask')

ngrid = np.array(maskgre)
ngrid=np.hstack(ngrid)
j=0
for i in range(len(ngrid)):
    if np.isfinite(ngrid[i]):
        ngrid[i]=j
        j=j+1
ngrid_gre=ngrid.reshape(maskgre.shape)

#%% GLaciers mask
# ds=xr.open_dataset('/Volumes/LaCie_NIOZ/data/barystatic/masks/barystatic_mask2.nc')
maskgla_rgi=np.array(ds.gla)
lat=np.array(ds.lat);lon=np.array(ds.lon)
dimlat=len(lat);dimlon=len(lon)
# unc.plot_world(lon,lat,maskgla_rgi,cmin=0,cmax=np.nanmax(maskgla_rgi),
#              cmap='tab20',title='Glaciers Mask')

maskgla = np.array(ds.mask3)
maskgla[maskgla<4]=np.nan
maskgla[np.isfinite(maskgla)]=1

maskgla=np.array(maskgla_rgi*maskgla)
maskgla[np.isfinite(maskgla)]=1
# unc.plot_world(lon,lat,maskgla,title='GLA mask')

ngrid = np.array(maskgla)
ngrid=np.hstack(ngrid)
j=0
for i in range(len(ngrid)):
    if np.isfinite(ngrid[i]):
        ngrid[i]=j
        j=j+1
ngrid_gla=ngrid.reshape(maskgla.shape)

#%% LWS mask
# ds=xr.open_dataset('/Volumes/LaCie_NIOZ/data/barystatic/masks/barystatic_mask2.nc')
mask_lws=np.array(ds.mask3)
mask_lws[np.where(mask_lws!=1)]=np.nan

lat=np.array(ds.lat);lon=np.array(ds.lon)
dimlat=len(lat);dimlon=len(lon)

# unc.plot_world(lon,lat,mask_lws,title='LWS mask')

ngrid = np.array(mask_lws)
ngrid=np.hstack(ngrid)
j=0
for i in range(len(ngrid)):
    if np.isfinite(ngrid[i]):
        ngrid[i]=j
        j=j+1
ngrid_lws=ngrid.reshape(mask_lws.shape)

#%% Open hector regional trends
t0=2003
t1=2016
path = '/Volumes/LaCie_NIOZ/data/barystatic/revisions/results_final/{}-{}/'.format(t0,t1)
# path = '/Volumes/LaCie_NIOZ/data/barystatic/results/{}-{}/'.format(t0,t1)
# ds=xr.open_dataset(path+'source_trend_temporal_unc_{}-{}.nc'.format(t0,t1))
ds=xr.open_dataset(path+'source_trend_unc_{}-{}.nc'.format(t0,t1))

print(ds)
names = list(np.array(ds.name))
name_ant = [name for name in names if name.startswith('AIS')]
name_gre = [name for name in names if name.startswith('GIS')]
name_lws = [name for name in names if name.startswith('LWS')]
name_gla = [name for name in names if name.startswith('GLA')]

lon=np.array(ds.lon)
lat=np.array(ds.lat)
dimlat=len(lat);dimlon=len(lon)
llon,llat=np.meshgrid(lon,lat)
llon=llon.flatten()
llat=llat.flatten()

grid_area=(sl.get_grid_area(np.zeros((mask.shape)))).reshape(180,360)
# unc.plot_world(lon,lat,grid_area,cmin=100,cmax=12000,clabel='Area (km2)',title='Grid Area')

df_scale= pd.DataFrame(llon,columns=['lon'])
df_scale['lat']=llat

df_slf = pd.DataFrame(llon,columns=['lon'])
df_slf['lat'] = - llat

df_tr = pd.DataFrame(llon,columns=['lon'])
df_tr['lat'] = - llat

#%% mask dict
mask_dict = {'AIS_IMB': maskant_regions,
             'AIS_UCI':maskant_basins,
             'AIS_CSR':maskant,
             'AIS_JPL':maskant,
             'GIS_CSR':maskgre,
             'GIS_JPL':maskgre,
             'GIS_IMB':maskgre_dyn,
             'GIS_UCI':maskgre_basins,
             'LWS_CSR':mask_lws,
             'LWS_JPL':mask_lws,
             'LWS_GWB':mask_lws,
             'LWS_WGP':mask_lws,
             'GLA_ZMP':maskgla,
             'GLA_WGP':maskgla,
             'GLA_CSR':maskgla,
             'GLA_JPL':maskgla,
             
             }

# import pickle
# def save_dict(obj, name, path ):
#     with open(path + name + '.pkl', 'wb') as f:
#         pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

# name = 'masks_dict_spatial_norm'
# pwd='/Volumes/LaCie_NIOZ/data/barystatic/revisions/results/'
# save_dict(mask_dict,name,pwd)
#%% # Distribute 1 mm and get spatial scale
for iname, name in enumerate(names):
    reg = name.split('_')[0]
    dataset = name.split('_')[1]
    mask=mask_dict[name]
    
    mask[np.isfinite(mask)] = 1
    da=ds.sel(name=[name])

    ic_idx =-1
    tr = np.array(da.best_trend[0,ic_idx,:,:])
    # unc.plot_ant(da.lon,da.lat,tr,cmin=-10,cmax=10,
    #          cmap=cm.cm.balance,
    #          title=str(name)+'\n trend',clabel='mm EWH/yr')

    # scale = unc.scale_pattern(tr,mask)
    scale = np.array(tr*mask).flatten()
    # if reg =='AIS':
    #     unc.plot_ant(lon,lat, scale.reshape(dimlat,dimlon),
    #              cmap=cm.cm.thermal,
    #              title=str(name)+'\nscaling pattern',cmin=-0.01,cmax=0.01,
    #             clabel=' scaling pattern')
    # if reg =='GIS':
    #     unc.plot_gre(lon,lat, scale.reshape(dimlat,dimlon),
    #              cmap=cm.cm.thermal,
    #              title=str(name)+'\nscaling pattern',cmin=-0.01,cmax=0.01,
    #             clabel=' scaling pattern')
    # else:
    #     unc.plot_world(lon,lat, scale.reshape(dimlat,dimlon),
    #              cmap=cm.cm.thermal,
    #              title=str(name)+'\nscaling pattern',cmin=-0.01,cmax=0.01,
    #             clabel=' scaling pattern')
    # print(unc.stats(scale))
    
    df_scale[name]=scale
#%% #%% run SLE
for iname, name in enumerate(names):
    slc = np.array(df_scale[name])
    # reg='gre'
    # transform in EWH:
    # sl_EWH = sle.height_to_EWH(slc)
    sl_EWH = np.array(slc)
    if name.split('_')[1]=='IMB':
        sl_EWH = -sl_EWH
    # input it in the SLE: 
    slf = sle.run_SLE(sl_EWH.reshape(dimlat,dimlon),str(name))
    # add it to df:
    df_slf[name] = slf.flatten()
    Z= np.array(df_slf[name]).reshape(dimlat,dimlon)
    # unc.plot_world(lon,-lat,Z,
    #                 cmin=-3,cmax=3,title=name,cmap=cmf.roma_r,
    #                 clabel='')

#%% open etopo to get ocean mask
ds_mask = xr.open_dataset('/Volumes/LaCie_NIOZ/data/ETOPO/ETOPO1_Ice-180x360.nc')
# ds_mask.z.plot(# vmin=-1,vmax=1
#                 );plt.show()

ds_mask=ds_mask.sortby('lat',ascending=False)
oceanmask=np.array(ds_mask.z)
oceanmask[oceanmask>=0]=np.nan
oceanmask[oceanmask<=0]=1
# plt.pcolor(oceanmask);plt.show()

# # oceanmask=oceanmask.flatten()
df2=df_slf[['lon','lat']]
for name in names:
    df2[name]=df_slf[name]*oceanmask.flatten()
# df2[names].hist()  

#%% Standard deviation of fingerprints
# ANT
df_ant = df_slf[name_ant]
std_ant=df_ant.std(axis=1)
#% %
# unc.plot_world(lon,-lat,np.array(std_ant).reshape(180,360),
#            cmin=0,cmax=2,
#            title='ANT \n Standard Deviation',
#            # fillcont=False,
#             cmap=cmf.devon_r,
#             # cmap=cmf.bilbao,          
#             # cmap=cmf.lajolla,
#            clabel='')
print(unc.stats(std_ant))

# GRE
df_gre = df_slf[name_gre]
std_gre=df_gre.std(axis=1)
#% %
# unc.plot_world(lon,-lat,np.array(std_gre).reshape(180,360),
#            cmin=0,cmax=2,
#            title='GRE \nStandard Deviation',
#            # fillcont=False,
#             cmap=cmf.devon_r,
#             # cmap=cmf.bilbao,          
#             # cmap=cmf.lajolla,
#            clabel='')
print(unc.stats(std_gre))

# GLA
df_gla = df_slf[name_gla]
std_gla=df_gla.std(axis=1)
#% %
# unc.plot_world(lon,-lat,np.array(std_gla).reshape(180,360),
#            cmin=0,cmax=2,
#            title='GLA \nStandard Deviation',
#            # fillcont=False,
#             cmap=cmf.devon_r,
#             # cmap=cmf.bilbao,          
#             # cmap=cmf.lajolla,
#            clabel='')
print(unc.stats(std_gla))

# LWS
df_lws = df_slf[name_lws]
std_lws=df_lws.std(axis=1)
#% %
# unc.plot_world(lon,-lat,np.array(std_lws).reshape(180,360),
#            cmin=0,cmax=2,
#            title='LWS \nStandard Deviation',
#            # fillcont=False,
#             cmap=cmf.devon_r,
#             # cmap=cmf.bilbao,          
#             # cmap=cmf.lajolla,
#            clabel='')
print(unc.stats(std_lws))



#%% make unc only normalized std
df_unc_norm=df_slf[['lon','lat']]
df_unc_norm['AIS']=std_ant
df_unc_norm['GIS']=std_gre
df_unc_norm['GLA']=std_gla
df_unc_norm['LWS']=std_lws

# save
pwd='/Volumes/LaCie_NIOZ/data/barystatic/revisions/results_final/'
df_unc_norm.to_pickle(pwd+"/not_normalized_spatial_unc.p")

# df_slf.to_pickle(pwd+"normalized_SLF.p")

#%%

# import xarray as xr
# import pandas as pd
# import numpy as np
# import string


# # plotting
# from cartopy import crs as ccrs , feature as cfeature
# import matplotlib.pyplot as plt
# # from matplotlib.gridspec import GridSpec
# from matplotlib.colors import ListedColormap
# # from matplotlib import colors
# import cmocean as cm
# # from cmcrameri import cm as cmf


# # user defined functions
# import sys
# sys.path.append("/Users/ccamargo/Documents/py_scripts/OM/")
# import utils_OM as om

# # warnings
# import warnings
# # With np.nanmean, for all-NaN slices, NaN is returned and a RuntimeWarning is raised.
# warnings.filterwarnings("ignore","Mean of empty slice", RuntimeWarning)
# warnings.filterwarnings("ignore","invalid value encountered in less", RuntimeWarning)
# warnings.filterwarnings("ignore","invalid value encountered in greater", RuntimeWarning)

# # %% define plotting preferences
# cmap_unc=cm.cm.curl
# # path_to_figures = '/Users/ccamargo/Documents/PhD/Barystatic/manuscript/revision_round1/figures/'

# dpi=300
# wi=20;hi=15
# dimlat=180;dimlon=360
# fontsize=25
# ticksize=20

# letters = list(string.ascii_lowercase)
# landcolor='darkgrey'

# # make our color scale
# col_dict={1:"black", # WN
#           2:"palegoldenrod", # PL
#           3:"lightpink", # PLWN
#           4:"orange", # AR1
#           5:"teal", # Ar5
#           6:"darkmagenta", # AR9
#           7:"skyblue", # ARf
#           8:"crimson" # GGM
#           }

# # We create a colormar from our list of colors
# cmapnm = ListedColormap([col_dict[x] for x in col_dict.keys()])

# #%% load data
# t0=2003;t1=2016

# path = '/Volumes/LaCie_NIOZ/data/barystatic/revisions/results/{}-{}/'.format(t0,t1)
# # path = '/Users/ccamargo/Desktop/'
# ds=xr.open_dataset(path+'final_dataset_{}-{}.nc'.format(t0,t1))
# lon=np.array(ds.lon)
# lat=np.array(ds.lat)
# X=np.array(lon);X[-1]=360
# Y=np.array(-lat)

# t0_a=1993
# t1_a=2017
# path_a = '/Volumes/LaCie_NIOZ/data/barystatic/revisions/results/{}-{}/'.format(t0_a,t1_a)
# # path = '/Users/ccamargo/Desktop/'
# ds_a=xr.open_dataset(path_a+'final_dataset_{}-{}.nc'.format(t0_a,t1_a))



# #% % mask
# mask=xr.open_dataset('/Volumes/LaCie_NIOZ/PhD/Barystatic/GRACE/mascons/JPL/global/LAND_MASK.CRI_360x180.nc')
# # mask.land_mask.plot()
# mask=mask.sortby('lat',ascending=False)
# oceanmask=np.array(mask.land_mask)
# oceanmask[oceanmask>0]=np.nan
# oceanmask[oceanmask==0]=1
# #%%
# clim=1
# cmap=cmap_unc
# c0=0
# cmin=0
# cmax=clim
# cmap = cm.tools.crop(cmap, cmin,cmax,c0)


# interval=0.1

# X[-1]=360
# X[0]=0
# # spatial unc is the same for all AIS, GIS, LWS and GLA (no matter the dataset)
# names=[name for name in np.array(ds.name) if name.split('_')[1]=='CSR']
# # names = ['AIS', 'GIS', 'GLA', 'LWS']
# fig = plt.figure(figsize=(15,10), facecolor='w',dpi=dpi)
# for iname, name in enumerate(names):
#     da=ds.sel(name=[name])
#     data = np.array(da.spatial_unc_norm[0,:,:]* oceanmask)
#     # data = np.array(df_unc_norm[name] * oceanmask.flatten()).reshape(180,360)
#     ax1 = plt.subplot(2,2,iname+1, projection=ccrs.Robinson())
    
#     ax1.coastlines(resolution='110m', zorder=3,color=landcolor) # zorder=3 makes sure that no other plots overlay the coastlines
#     ax1.add_feature(cfeature.LAND,color=landcolor,# alpha=0.5,
#                     zorder=3)
#     ax1.set_global() # make sure the projection is maximised inside the plot. Makes a circle
    
#     lv=np.arange(cmin,cmax+interval,interval)
#     csf=plt.contourf(X,Y,np.abs(data),levels=lv,
#               transform = ccrs.PlateCarree(),cmap=cmap)
#     mu,glb,mask=om.reg_to_glb(data,Y,X)
#     glb=np.round(glb,2)
 
#     cs=ax1.contour(X,Y,data,
#                    levels=[glb],
#                # levels=[0.15, 0.30, 0.45, 0.6, 0.75],
#             # vmin=-0.6,
#             # vmax=0.6,
#         transform = ccrs.PlateCarree(),
#         #cmap='coolwarm',#extend='both'
#         colors=('black',),linestyles=('--',),linewidths=(2,)
#         )   
#     ax1.clabel(cs,cs.levels,fmt='%5.2f',colors='k',fontsize=12)
    
#     cp=plt.pcolormesh(X,Y,data,
#                 vmin=cmin,vmax=cmax,
#                 zorder=0,
#                 transform = ccrs.PlateCarree(),cmap=cmap)
#     ax1.set_title('({}). '.format(letters[iname])+str(name).split('_')[0]+'\n{}'.format(np.round(glb,3)),size=fontsize)    

    
# plt.tight_layout()
# cbar_ax2 = fig.add_axes([0.025, -0.05, 0.95, 0.04])
# cbar2 = plt.colorbar(csf, cax=cbar_ax2, orientation='horizontal')
# cbar2.set_label(label='Spatial-structural Uncertainty (mm/yr)\n  Normalized',size=fontsize, family='serif')
# cbar2.ax.tick_params(labelsize=ticksize) 

# plt.show()
    
# #%%
# # names=[name for name in np.array(ds.name) if name.split('_')[1]=='CSR']
# names = ['AIS', 'GIS', 'GLA', 'LWS']
# fig = plt.figure(figsize=(15,10), facecolor='w',dpi=dpi)
# for iname, name in enumerate(names):
#     # da=ds.sel(name=[name])
#     # data = np.array(da.spatial_unc_norm[0,:,:]* oceanmask)
#     data = np.array(df_unc_norm[name] * oceanmask.flatten()).reshape(180,360)
#     ax1 = plt.subplot(2,2,iname+1, projection=ccrs.Robinson())
    
#     ax1.coastlines(resolution='110m', zorder=3,color=landcolor) # zorder=3 makes sure that no other plots overlay the coastlines
#     ax1.add_feature(cfeature.LAND,color=landcolor,# alpha=0.5,
#                     zorder=3)
#     ax1.set_global() # make sure the projection is maximised inside the plot. Makes a circle
    
#     lv=np.arange(cmin,cmax+interval,interval)
#     csf=plt.contourf(X,Y,np.abs(data),levels=lv,
#               transform = ccrs.PlateCarree(),cmap=cmap)
#     mu,glb,mask=om.reg_to_glb(data,Y,X)
#     glb=np.round(glb,2)
 
#     cs=ax1.contour(X,Y,data,
#                    levels=[glb],
#                # levels=[0.15, 0.30, 0.45, 0.6, 0.75],
#             # vmin=-0.6,
#             # vmax=0.6,
#         transform = ccrs.PlateCarree(),
#         #cmap='coolwarm',#extend='both'
#         colors=('black',),linestyles=('--',),linewidths=(2,)
#         )   
#     ax1.clabel(cs,cs.levels,fmt='%5.2f',colors='k',fontsize=12)
    
#     cp=plt.pcolormesh(X,Y,data,
#                 vmin=cmin,vmax=cmax,
#                 zorder=0,
#                 transform = ccrs.PlateCarree(),cmap=cmap)
#     ax1.set_title('({}). '.format(letters[iname])+str(name).split('_')[0]+'\n{}'.format(np.round(glb,3)),size=fontsize)    

    
# plt.tight_layout()
# cbar_ax2 = fig.add_axes([0.025, -0.05, 0.95, 0.04])
# cbar2 = plt.colorbar(csf, cax=cbar_ax2, orientation='horizontal')
# cbar2.set_label(label='Spatial-structural Uncertainty (mm/yr)\n  Not Normalized',size=fontsize, family='serif')
# cbar2.ax.tick_params(labelsize=ticksize) 

# plt.show()
    
