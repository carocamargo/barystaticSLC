#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  9 08:33:45 2021

Functions used for analysis and figures of the 
Ocean mass (barystatic) regional sea-level change manuscript

@author: ccamargo
"""


import pandas as pd
import xarray as xr
import numpy as np
import pickle

# plotting
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


#%% user defined functions:
def load_dict(name, path ):
    with open(path + name + '.pkl', 'rb') as f:
        return pickle.load(f)
    
#%%   
def quadrsum(X):
    if len(X.shape)==3:
        Y = np.zeros((X.shape[1],X.shape[2]))
    else:
        Y = np.zeros((X.shape[1]))
    for i in range(X.shape[0]):
        Y = Y + (X[i]*X[i])
    Y = np.sqrt(Y)
    return Y
#%%
def get_percentage_quadra(X):
    y_sum = quadrsum(X) * quadrsum(X) # remove the square root
    x_rel = (X*X)/y_sum    
    return x_rel
#%%
def get_percentage(X):
    x_rel = X/np.nansum(X,axis=0)
    return x_rel
#%%

def plot_pie_inset(data,ilon,ilat,ax,width,
                   colors=['seagreen','purple','cyan'],
                  ):
    ax_sub= inset_axes(ax, width=width, height=width, loc=10, 
                       bbox_to_anchor=(ilon, ilat),
                       bbox_transform=ax.transData, 
                       borderpad=0)
    wedges,texts= ax_sub.pie(data,colors=colors,counterclock=False,
                                 startangle=90)

    ax_sub.set_aspect("equal")
#%%
def plot_clustered_stacked(dfall, labels=None,
                           dpi=300,
                           fsize=(15,10),
                           fontsize=15,
                           ylabel='',
                           xoffset=0,
                           sideoffset=0,
                           title="multiple stacked bar plot",  
                           xlabel=False,
                           H="/", **kwargs):
    """
    Given a list of dataframes, with identical columns and index, 
    create a clustered stacked bar plot. 
    labels is a list of the names of the dataframe, used for the legend
    title is a string for the title of the plot
    H is the hatch used for identification of the different dataframe
    """

    n_df = len(dfall) # number of dataframes
    n_col = len(dfall[0].columns) #number of variables in each df
    n_ind = len(dfall[0].index) # x locations
    fig = plt.figure(dpi=dpi,
                     figsize=fsize)
    axe = plt.subplot(111)
    hatches = ['','**', 'o',
               '--',  'xx', 'oo', 'OO',  '**','\\', '||',]
    hatches = hatches[0:n_df] *n_df

    for df in dfall : # for each data frame

        axe = df.plot(kind="bar",
                      linewidth=1,
                      stacked=True,
                      ax=axe,
                      color=['seagreen','purple','cyan'],
                      legend=False,
                      grid=False,
                      #**kwargs
                      )  # make bar plots
        # plt.show()

    h,l = axe.get_legend_handles_labels() # get the handles we want to modify
    # h is each x-location
    # l is each variable (repeating for each df)
    for i in range(0, len(l), n_col): # for each variable
        for j, pa in enumerate(h[i:i+n_col]):
            for ir,rect in enumerate(pa.patches): # for each index
                rect.set_x(rect.get_x() + 1 / float(n_df + 1) * i / float(n_col))
                # rect.set_hatch(hatches[i]) #hatch according to each df
                # for jj in np.arange(n_df):
                #     rect.set_hatch(hatches[jj])
                rect.set_width(1 / float(n_df + 1))
    for i in range(len(h[0])):
        rect=h[0][i]
        for ilb, lb in enumerate(labels):
            # height = rect.get_height()
            height = xoffset
            x= rect.get_x()+0.05 +(ilb/5) + sideoffset
            plt.text(x, 
                     height,

                str(lb),fontsize=15,
                rotation = 90,
                ha='left', va='bottom')
    #% %
    axe.set_xticks((np.arange(0, 2 * n_ind, 2) + 1 / float(n_df + 1)) / 2.)
    if xlabel:
        axe.set_xticklabels(xlabel, rotation = 25,fontsize=fontsize)

    else:
        axe.set_xticklabels(df.index, rotation = 0,fontsize=fontsize)
    axe.set_title(title,fontsize=fontsize+5)
    axe.tick_params(axis='both', which='major', labelsize=fontsize)
    axe.tick_params(axis='x', which='major', pad=30)
    axe.set_ylabel(ylabel,fontsize=fontsize)
    # axe.set_xlabel(xlabel,labelpad=10)
    axe.set_xlim([-0.5,n_ind-0.25])
    # Add invisible data to add another legend
    n=[]        
    for i in range(n_df):
        n.append(axe.bar(0, 0, color="lightgray", 
                         hatch= hatches[i]# H * i
                         ))

    l1 = axe.legend(h[:n_col], l[:n_col], loc=[1.01, 0.5],fontsize=fontsize)
    # Add labels above the each bar graphs



    # if labels is not None:
    #     l2 = plt.legend(n, labels, loc=[1.01, 0.1],fontsize=fontsize) 
    axe.add_artist(l1)
    # return fig

    return fig


#%%
def from_360_to_180(lon_in):
    # given a longitude array that goes from 0 to 360, 
    # returns an array that goes from -180 to 180
    lon_out=np.copy(lon_in)
    for i,ilon in enumerate(lon_in):
        if ilon > 180: 
            lon_out[i]=ilon-360
    return lon_out


    #%%
def get_ocean_area(lat,lon,dry_mask):
    # Function that computes the ocean area given a latitude and longitude arrays.
    # Latitude and longitude should be in decimal degrees
    # It requires a dry mask: an array with 0 for where is land and 1 where is ocean.
    # Returns the ocean surface area in kmˆ2 (float number)
    # It also returns an area matrix (km^2), of size(len(lat),len(lon)). 
    # This is the 'weight' of the area in each grid cell. 
    
    R=6371 # Earth's radius in km
    
    #Check the grid resolution:
    deltalat=180/len(lat);
    deltalon=360/len(lon) 
    
    #Transform from degrees to km:

    deltay=(2*np.pi*R*deltalat)/360 #lat to km
    deltax=(2*np.pi*R*np.cos(np.radians(lat))*deltalon)/360 #lon to km
    
    area=np.array([deltax*deltay]*len(lon)).transpose()
    ocean_surf=np.sum(area*dry_mask)

    
    return ocean_surf, area
#%%
def reg_to_glb(value,lat,lon,*args,
               save=False,path='',buf=False,
               path_mask='/Users/ccamargo/Documents/PhD/Barystatic/SLM_run/model/topo/'):
    #value is a 4d (depth,time,lat,lon), 3d (time,lat,lon) or a 2D(lat,lon) matrix
    # ATT: lat and lon should be the last positional arguments of the matrix!!
    # lat and lon are mandatory, because we have a regional grid
    # time and depth are optional arguments that should be provided in case or a 3D or 4D matrix
    #dimensions of the output will depend on the dimensions of the input value matrix
    # fmask= D : computes land mask based on the datset (assumes that the dataset has nan over land)
    # fmask=JPL-l : uses land mask of JPL'/Users/ccamargo/Documents/PhD/Barystatic/GRACE/mascons/JPL/global/LAND_MASK.CRI_360x180.nc'
    #         -l: uses 0=land,1=ocean (return changes over the ocean)
    #         -o: uses 1=land, 0=ocean (return changes over land)
    #================
    # Subroutines
    #===============   
    #######################################
    # Check dimensions of value matrix
    s=value.shape 
    
    #make an array to compare with shape of value matrix
    d=np.array([len(lat),len(lon)])
    if args:
        for i,n in enumerate(args):
            d=np.append(d,len(n))
            
    # check if dimensions match:
    if any(s)!=any(d):
        raise ValueError('Error: dimensions do not agree')
       
    ####################################### 
    # Get ocean-land mask:
    dimlat=np.array(lat.shape)[0]
    dimlon=np.array(lon.shape)[0]
    f=path_mask+'/mask-'+str(dimlat)+'.xyz'
    df=pd.read_csv(f,header=None,sep='\s+')
    df.columns=['lon','lat','mask']
    oceanmask=np.array(df['mask']).reshape(dimlat,dimlon) # 1 for land and 0 for ocean

    ####################################### 
    
    # Get ocean area based on mask:
    surf,area = get_ocean_area(lat,lon,oceanmask)
    surf_m=surf*1000000 # m2
    area_m=area*1000000 # m2

    #=========================
    # Main programme
    #=========================   
    # To go from regional to global we need first:
    # multiply by the area of each grid, because grid cells do not have the same size
    # and multiply by the land mask, which will garantee value 0 in the land grids
    # Then divide by the entire ocean surface area:

    #glb=np.sum(value*area_m*mask)/surf_m
    da=xr.Dataset(data_vars={'data':(('lat','lon'),value),
                             'oceanmask':(('lat','lon'),oceanmask)
                             },
                     coords={'lat':lat,
                             'lon':lon})
    mu=(da.data*area_m*oceanmask).sum(dim=('lat','lon'))/surf_m
    glb=np.array(mu.data)

    return mu, glb, oceanmask

def get_grid_area(grid):
    # given a grid that has dimensions: len(lat) x len(lon), compute the 
    #area in each grid cell, in km2
    # input: grid: numpy array
    
    earth_radius = 6371 # km
    # earth_radius = 6378137/1000# a more precise earth radius
    earth_diam = 2* earth_radius # diameter in km
    earth_circ = np.pi* earth_diam # earth's circunference in meters
    
    #% %
    dimlat=grid.shape[0]
    dimlon=grid.shape[1]
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



