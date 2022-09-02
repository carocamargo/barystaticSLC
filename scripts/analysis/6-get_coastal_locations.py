#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  7 17:54:18 2022

@author: ccamargo
"""


import pandas as pd
import xarray as xr
from cartopy import crs as ccrs# , feature as cfeature
import matplotlib.pyplot as plt
import numpy as np

#%% find closest points of our TG, make sure that they are water
def harvedist(lats,lons,qlats,qlons):
    # https://www.movable-type.co.uk/scripts/latlong.html
    #calculate angular distances between query coordinates (qlats,qlons) and data array grid (lats,lons)
    lat0,lat = np.meshgrid(np.radians(lats),np.radians(qlats))
    lon0,lon = np.meshgrid(np.radians(lons),np.radians(qlons))
    # lat0,lat = np.meshgrid(lats,qlats)
    # lon0,lon = np.meshgrid(lons,qlons)
    
    lat=np.radians(lat);lat0=np.radians(lat0)
    lon=np.radians(lon);lon0=np.radians(lon0)
    delta_lat = np.array(lat-lat0)
    delta_lon = np.array(lon-lon0)
    R =6373 # km
    a = np.sin(delta_lat / 2)**2 + np.cos(lat) * np.cos(lat0) * np.sin(delta_lon / 2)**2

    c = 2 * np.arctan2(np.sqrt(a),np.sqrt(1-a))
    
    # d = R * c
    # print(d)
    return(np.degrees(c))

# def find_nearest(lats,lons,qlats,qlons):
#     #finds nearest coordinates to query latitude
#     #and longitude pairs based on minimum angular distance [deg]
        
#     #calculate angular distances between query coordinates and data array grid
#     dists,dist_km =harvedist(lats,lons,qlats,qlons)
    
#     min_dists = np.nanmin(dists,axis=1) #find minimum angular distances in grid to query points
#     min_idx = np.nanargmin(dists,axis=1) #indices
#     # min_dist_km = dist_km[min_idx]
    
#     out_lat = lats.flatten()[min_idx]
#     out_lon = lons.flatten()[min_idx]
    
#     return min_dists, min_idx, out_lat, out_lon
#% %
def angdist(lats,lons,qlats,qlons):
    lat0,lat = np.meshgrid(np.radians(lats),np.radians(qlats))
    lon0,lon = np.meshgrid(np.radians(lons),np.radians(qlons))
    
    temp = np.arctan2(np.sqrt((np.cos(lat)*np.sin(lon-lon0))**2 + (np.cos(lat0)*np.sin(lat) - np.sin(lat0)*np.cos(lat) * np.cos(lon-lon0))**2),
                      (np.sin(lat0)*np.sin(lat) + np.cos(lat0)*np.cos(lat)*np.cos(lon-lon0)))
    return(np.degrees(temp))

def find_nearest(da,qlats,qlons):
    #finds nearest unmasked ocean grid cell in xarray dataarray to query latitude
    #and longitude pairs based on minimum angular distance [deg]
    
    #fetch coordinate names
    lonname = np.array(da.coords)[['lon' in x for x in da.coords]][0]
    latname = np.array(da.coords)[['lat' in x for x in da.coords]][0]

    #get lats & lons
    lats = np.array(da[latname])
    lons = np.array(da[lonname])

    if lats.shape!=lons.shape:
        lats,lons = np.meshgrid(lats,lons)
   
        
    #calculate angular distances between query coordinates and data array grid
    dists=angdist(lats,lons,qlats,qlons)
    
    #mask land out
    if 'time' in da.dims:
        dists[0,~np.isfinite(da.isel(time=0).values.flatten())] = np.nan
    else:
        dists[0,~np.isfinite(da.values.flatten())] = np.nan
    
    min_dists = np.nanmin(dists,axis=1) #find minimum angular distances in grid to query points
    min_idx = np.nanargmin(dists,axis=1) #indices
    
    
    out_lat = lats.flatten()[min_idx]
    out_lon = lons.flatten()[min_idx]
    #potentially build in a filter here if unreasonably large distances

    return min_dists, min_idx, out_lat, out_lon



#%%
periods =[
    # (2005,2016),
          (1993,2018),
          (2003,2017)
          ]
for period in periods:
    t0=period[0]
    t1=period[1] -1
    #% % open OM dataset
    # t0=2005
    # t1=2015
    path = '/Volumes/LaCie_NIOZ/data/barystatic/revisions/results_final/{}-{}/'.format(t0,t1)
    # path = '/Users/ccamargo/Desktop/'
    df=pd.read_pickle(path+'OM_reconstructions_{}-{}.p'.format(t0,t1))
    
    
    
    #% % sel 10
    cat = pd.read_pickle('/Volumes/LaCie_NIOZ/data/barystatic/coastal_locations.p')
    
    cities= [
        'Rotterdam', 
        'Rio_de_Janeiro',
         'Jakarta',
             'Tokyo',
             'Lima',
             'Cape Town',
              'Muqdisho (Mogadishu)', # Somalia
            'Sydney',
            'Vancouver',
            
            'Virginia Beach', 
            
            ]
    cat = cat[cat['City'].isin(cities)].reset_index()
    
    
    #% % add mask
    ds=xr.open_dataset('/Volumes/LaCie_NIOZ/data/barystatic/masks/LAND_MASK_CRI-JPL_180x360_conservative.nc')
    ds.mask.plot()
    ds=ds.sortby('lat',ascending=False)
    # % % add mask
    # path_mask='/Users/ccamargo/Documents/PhD/Barystatic/SLM_run/model/topo/'
    # dimlat=180;dimlon=360
    # f=path_mask+'/mask-'+str(dimlat)+'.xyz'
    # df_mask=pd.read_csv(f,header=None,sep='\s+')
    # df_mask.columns=['lon','lat','msk']
    # df['msk']=df_mask['msk']
    #% %
    # df['msk']=np.array(ds.mask).flatten()
    mask=np.array(ds.mask)
    mask[130:150,4:5]=10;
    print(ds.lat[130:150])
    df['msk']=mask.flatten()
    df['llat']=df['lat']
    df['llon']=df['lon']
    df_multiindex = df.set_index([ 'lat','lon'])
    ds=df_multiindex.to_xarray()
    
    
    # ds.msk.plot()
    # fig = plt.figure(dpi=300)
    # ax = plt.axes(projection=ccrs.PlateCarree())
    # plt.title('mask')
    # ax.coastlines(resolution='110m')
    # ax.set_global()
    # # plt.savefig('map.png')
    # splot = ax.scatter(
    #     "lon","lat",
    #     c="msk", # color by a variable
    #     data=df,
    #     s=10,
    #     cmap="plasma",
    #     transform=ccrs.PlateCarree(),
    #     zorder=0
    #     )
    # # plt.colorbar(splot, orientation ='horizontal', label='OM Trend')
    # plt.show()
    
    
    #% % find closest point in the 1degree resolution
    da=ds.msk
    min_dist,min_idx,out_lat,out_lon= find_nearest(da,
                                                   cat['lat'],cat['lon']
                                                   )
    #% %
    # cat['lat_new']=out_lat
    # cat['lon_new']=out_lon
    stations=pd.DataFrame(out_lon,columns=['lon'])
    stations['lat']=out_lat
    stations['City']=cat['City']
    stations['Country']=cat['Country']
    
    # merge with mask and remove any station located on land
    df2=pd.merge(stations,df,on=['lon','lat'])
    # df2.drop(df2[df2.msk<1].index,inplace=True)
        
    #% % plot 
    # fig = plt.figure(dpi=300)
    # ax = plt.axes(projection=ccrs.PlateCarree())
    # plt.title('Coastal stations')
    # ax.coastlines(resolution='110m')
    # ax.set_global()
    # # plt.savefig('map.png')
    # splot = ax.scatter(
    #     "lon","lat",
    #     c="msk", # color by a variable
    #     data=df2,
    #     s=20,
    #     cmap="plasma",
    #     transform=ccrs.PlateCarree(),
    #     zorder=0
    #     )
    # plt.colorbar(splot, orientation ='horizontal', label='OM Trend')
    # plt.show()
    #% % save df
    df2.to_pickle(path+'coastal_examples_10_{}-{}.p'.format(t0,t1))
