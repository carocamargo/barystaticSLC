#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  6 15:00:00 2022

@author: ccamargo
"""

import xarray as xr
import pandas as pd
import numpy as np
import string


# plotting
from cartopy import crs as ccrs , feature as cfeature
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import ListedColormap
from matplotlib import colors
import cmocean as cm
from cmcrameri import cm as cmf


# user defined functions
import sys
sys.path.append("/Users/ccamargo/Documents/py_scripts/OM/")
import utils_OM as om

# warnings
import warnings
# With np.nanmean, for all-NaN slices, NaN is returned and a RuntimeWarning is raised.
warnings.filterwarnings("ignore","Mean of empty slice", RuntimeWarning)
warnings.filterwarnings("ignore","invalid value encountered in less", RuntimeWarning)
warnings.filterwarnings("ignore","invalid value encountered in greater", RuntimeWarning)
#%% define plotting preferences
cmap_unc=cm.cm.curl
path_to_figures = '/Users/ccamargo/Documents/PhD/Barystatic/manuscript/revision_round1/revision-submission/figures/'

dpi=300
wi=20;hi=15
dimlat=180;dimlon=360
fontsize=25
ticksize=20

letters = list(string.ascii_lowercase)
landcolor='darkgrey'

# make our color scale
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

#%% load data
t0=2003;t1=2016
path = '/Volumes/LaCie_NIOZ/data/barystatic/results/{}-{}/'.format(t0,t1)
ds=xr.open_dataset(path+'final_dataset_{}-{}.nc'.format(t0,t1))

path = '/Volumes/LaCie_NIOZ/data/barystatic/revisions/results_final/{}-{}/'.format(t0,t1)
# path = '/Users/ccamargo/Desktop/'
ds2=xr.open_dataset(path+'final_dataset_{}-{}.nc'.format(t0,t1))

lon=np.array(ds.lon)
lat=np.array(ds.lat)
X=np.array(lon);X[-1]=360
Y=np.array(-lat)

#% % mask
mask=xr.open_dataset('/Volumes/LaCie_NIOZ/PhD/Barystatic/GRACE/mascons/JPL/global/LAND_MASK.CRI_360x180.nc')
# mask.land_mask.plot()
mask=mask.sortby('lat',ascending=False)
oceanmask=np.array(mask.land_mask)
oceanmask[oceanmask>0]=np.nan
oceanmask[oceanmask==0]=1
#%%
nrows=4
ncols=4
fig = plt.figure(figsize=(15,10), facecolor='w')
interval=0.1

gs = GridSpec(nrows, ncols)
irow = 0

alpha=0.9

regions = ['GLA','LWS','GLA','LWS']
combos= ['JPL','CSR','JPL','CSR']
panels=['Trend', 'Uncertainty']
name_list = [
             'GLA_JPL','LWS_JPL',
             'GLA_CSR','LWS_CSR',
             'GLA_JPL','LWS_JPL',
             'GLA_CSR','LWS_CSR',
             
             ]
# ncols=len(panels)
iletter=0
for iname, name in enumerate(name_list):
    if iname<4:
        
        da=ds.sel(name=name)
    else:
        da=ds2.sel(name=name)
    
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
            Y=np.array(-lat);X=np.array(lon);X[0]=0;X[-1]=360
            clim=1
            cmap=cmf.roma_r
            data = np.array(da.trend[:,:])  


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
            data = np.abs(da.uncs[0,:,:])

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

kurs=path_to_figures+'mascon_split.png'
        
fig.savefig(kurs,format='png',dpi=300,bbox_inches='tight')
