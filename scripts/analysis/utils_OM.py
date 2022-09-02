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
import scipy.stats as st
from alive_progress import alive_bar

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

#%%
def design_matrix(time,acc_idx):
    #Build Design Matrix:
    A=np.ones([np.sum(acc_idx),7])
    A[:,1]=time[acc_idx]-np.nanmean(time) #first degree term
    A[:,2] = 0.5 * (time[acc_idx] - np.mean(time[acc_idx]))**2 #second degree term
    A[:,3]=np.cos(2*np.pi*time[acc_idx]) # cos-annual
    A[:,4]=np.sin(2*np.pi*time[acc_idx]) # sin-annual
    A[:,5]=np.cos(2*np.pi*2*time[acc_idx]) #cos-semi-annual
    A[:,6]=np.sin(2*np.pi*2*time[acc_idx]) #sin semi-annual
    return A

def OLS(x,y):
    time=x
    height=y
    acc_idx = np.isfinite(height)
    
    ##desgin matrix
    A  =design_matrix(time, acc_idx)
    #Fit OLS to the matrix
    solution=np.linalg.lstsq(A,y[acc_idx],rcond=None)[0]
    trend=solution[1] # trend

    return trend

def create_noise(mu,sigma,t):
    noise = np.random.normal(mu,sigma,t)
    return noise

def max_range(data):
    return max(np.abs(np.nanmax(data) - np.nanmean(data)),
                        np.abs(np.nanmin(data) - np.nanmean(data)))

# def MC_perturb(time,rate,unc,N):
#     # Given a time series of an estimate and uncertainties
#     # Perturb the time series N times, with random noise*unc
#     # get the maximum distance from the mean as the intrinsic unc
#     acc_idx = np.isfinite(rate)
#     t = len(rate[acc_idx])
#     # samples = np.zeros((N,t))
#     slopes = np.zeros((N))
#     np.random.seed(1) # Fix random number generator

#     for i in range(N):
#         # noise = np.array(create_noise(0,1,t)*unc)
#         # samples[i] = np.array(rate[acc_idx] + create_noise(0,1,t)*unc[acc_idx])
#         slopes[i] = OLS(time[acc_idx], np.array(rate[acc_idx] + create_noise(0,1,t)*unc[acc_idx]))
        
#     # perturbed = np.nanmean(slopes)
#     perturbed_error = max_range(slopes)
#     return perturbed_error

def MC_perturb(time,rate,unc,N,
               ci_level=0.99,
               plot=False,
               unit='mm',
               dataset='',
               fontsize=15):
    # Given a time series of an estimate and uncertainties
    # Perturb the time series N times, with random noise*unc
    # get the maximum distance from the mean as the intrinsic unc
    acc_idx = np.isfinite(rate)
    t = len(rate[acc_idx])
    samples = np.zeros((N,t))
    slopes = np.zeros((N))
    np.random.seed(1) # Fix random number generator
    if plot:
        plt.figure(figsize=(10,15))
        ncol=2
        nrow=1
        plt.subplot(ncol,nrow,1)
    for i in range(N):
        # noise = np.array(create_noise(0,1,t)*unc)
        samples[i] = np.array(rate[acc_idx] + create_noise(0,1,t)*unc[acc_idx])
        slopes[i] = OLS(time[acc_idx], samples[i])
        if plot:
            plt.plot(time[acc_idx],samples[i],alpha=0.1,color='gray')
    # perturbed = np.nanmean(slopes)
    # perturbed_error = max_range(slopes)
    x=slopes
    ci = st.norm.interval(alpha=ci_level, loc=np.mean(x), scale=x.std())
    perturbed_error = max(np.abs(ci-slopes.mean()))
    ci_width = np.abs(ci[0]-ci[1])
    perturbed_error = ci_width/2
    
    if plot:
        plt.plot(time[acc_idx],rate[acc_idx],linewidth=2,color='black',label='rate')
        plt.plot(time[acc_idx],rate[acc_idx]+unc[acc_idx],color='red',alpha=0.5,label='rate ± unc')
        plt.plot(time[acc_idx],rate[acc_idx]-unc[acc_idx],color='red',alpha=0.5)
        plt.ylabel(unit,fontsize=fontsize)
        plt.xlabel('time',fontsize=fontsize)
        plt.legend(fontsize=fontsize)
        plt.title('{}: {} ± {} {}/year'.format(dataset,
                                            np.round(x.mean(),3),
                                            np.round(perturbed_error,3),
                                            unit),fontsize=fontsize)
        plt.subplot(ncol,nrow,2)
        #% %
        plt.hist(slopes, bins=20, color='c', edgecolor='k',
                          density=True,
                          alpha=0.65)
        plt.xlabel(unit,fontsize=fontsize)
        plt.axvline(slopes.mean(),color='k',linestyle='--')
        min_ylim, max_ylim = plt.ylim()
        
        plt.text(x.mean()*1.01, max_ylim*0.9, 'Mean: {:.3f}'.format(x.mean()),fontsize=fontsize)
        plt.axvline(ci[0],c='red',linestyle='--',alpha=0.5,label='{}% CI'.format(ci_level*100))
        plt.axvline(ci[1],c='red',linestyle='--',alpha=0.5)
        plt.legend(loc='upper left',fontsize=fontsize)
        plt.ylabel('Density distribution',fontsize=fontsize)
        
        plt.tight_layout()
        plt.show()
    return perturbed_error

def intrinsic(time,rate,unc,N,dim='3D',ci_level=0.95):
    
    # given a rate and unc arrays, that have dimensions
    # time,lat,lon,
    # compute the intrinsic uncertainty rate 
    # return array that has dimensions latxlon
    if dim=='3D':
        dimtime,dimlat,dimlon = rate.shape
        rate = rate.reshape(dimtime,dimlat*dimlon)
        unc = unc.reshape(dimtime,dimlat*dimlon)
        error = np.zeros((dimlat*dimlon))
        error.fill(np.nan)
        j=0
        # start = dt.now()
        with alive_bar(dimlat*dimlon) as bar:
            for i in range(dimlat*dimlon):
                if np.any(np.isfinite(rate[:,i])):
                    j=j+1
                    # start = dt.now()
                    error[i] = MC_perturb(time, rate[:,i], unc[:,i], N)
                bar()
                    # deltaT = dt.now()-start
                    # print(deltaT)
        error = error.reshape(dimlat,dimlon)
    else:
        error = MC_perturb(time,rate,unc,N,ci_level=ci_level)

    return error

#%%
def gaus_filter(U, sigma=0.8, truncate=4.0):
    """
    Applies a gaussian filter to a numpy array that is not disturbed by NaNs.
    The code is adapted from an answer to this question:
    https://stackoverflow.com/questions/18697532/gaussian-filtering-a-image-with-nan-in-python
    This code uses the gaussian filter from scipy image. It calculates the 
    kernel size internally based on sigma and the truncate parameters as
    int(truncate * sigma + 0.5).

    Parameters
    ----------
    U : numpy array
        Array of the data to which the filter shall be applied.
    sigma : float, optional
        Standard deviation of the gaussian filter. The default is 0.8.
    truncate : float, optional
        Truncate filter at this many sigmas. The default is 4.0.

    Returns
    -------
    Z : TYPE
        DESCRIPTION.

    """
    import numpy as np
    from scipy.ndimage import gaussian_filter
       
    V = U.copy()
    V[np.isnan(U)] = 0
    VV = gaussian_filter(V,sigma=sigma,truncate=truncate)
    
    W = 0*U.copy()+1
    W[np.isnan(U)] = 0
    WW = gaussian_filter(W,sigma=sigma,truncate=truncate)
    
    # replace land with nan again to avoid invalid value in true divide
    WW[np.isnan(U)] = np.nan
    Z = VV/WW

    return Z    


