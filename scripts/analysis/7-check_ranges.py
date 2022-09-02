#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 12:16:09 2022

@author: ccamargo
"""

import pandas as pd
import numpy as np
import scipy.stats as st
import pickle
def load_dict(name, path ):
    with open(path + name + '.pkl', 'rb') as f:
        return pickle.load(f)
def open_mask():
    name = 'mask'
    path = '/Volumes/LaCie_NIOZ/data/barystatic/revisions/'
    
    dic = load_dict(name,path)
    mask = dic['mask']
    mask[mask!=0] = np.nan
    mask[np.isfinite(mask)] = 1
    return mask
# import numpy as np, scipy.stats as st

# returns confidence interval of mean
def confIntMean(a, conf=0.95):
  mean, sem, m = np.mean(a), st.sem(a), st.t.ppf((1+conf)/2., len(a)-1)
  return mean - m*sem, mean + m*sem
#%%
def stats(combo):
   # print('')
    print(combo)
    mask = open_mask()
    mask=mask.flatten()
    for var in ['trend','unc']:
        print(var)
        p5 = []
        p95 = []
        x = np.array(df['{}_{}_tot'.format(combo,var)]*mask)
        x = x[np.isfinite(x)]
        p5.append(np.round(np.percentile(x, 5),2))
        p95.append(np.round(np.percentile(x, 95),2))
        print('95%: {} to {}'.format(np.min(p5),np.max(p95)))
    return '\n'
   
#%% 95thrange
periods  = [(2003,2016),
             (1993,2017)
            ]
for period in periods:
    t0,t1 = period
    print('\nFrom {} to {}'.format(t0,t1))
    path='/Volumes/LaCie_NIOZ/data/barystatic/revisions/results_final/{}-{}/'.format(t0,t1)

    file = 'OM_reconstructions_{}-{}.p'.format(t0,t1)

    df = pd.read_pickle(path+file)

    if t0>2002:
        combos = ['JPL','CSR',
                      'IMB+WGP',
                       'IMB+GWB+ZMP',
                       'UCI+WGP','UCI+GWB+ZMP',
                      ]
    else:
        combos = [#'JPL','CSR',
                      'IMB+WGP',
                       'IMB+GWB+ZMP',
                       'UCI+WGP','UCI+GWB+ZMP',
                      ]   

    means = []
    mins = []
    maxs = []
    p5 = []
    p95 = []
    for combo in combos:
        # print('')
        # print(combo)
        # print('trend')
        mask = open_mask()
        mask=mask.flatten()
        x = np.array(df['{}_trend_tot'.format(combo)]*mask)
    
        x = x[np.isfinite(x)]
        ci = st.norm.interval(alpha=0.95, loc=np.mean(x), scale=x.std())
        means.append(np.round(np.nanmean(x),2))
        mins.append(np.round(np.nanmin(x),2))
        maxs.append(np.round(np.nanmax(x),2))
        p5.append(np.round(np.percentile(x, 5),2))
        p95.append(np.round(np.percentile(x, 95),2))
        
    print('95th percentile for all datasets:')
    print('trend:')
    print('95%: {} to {}'.format(np.min(p5),np.max(p95)))
    p5 = []
    p95 = []
    for combo in [#'JPL','CSR',
                  'IMB+WGP','IMB+GWB+ZMP',
                  'UCI+WGP','UCI+GWB+ZMP',]:
        # print('')
        # print(combo)
        # print('trend')
        mask = open_mask()
        mask=mask.flatten()
        x = np.array(df['{}_trend_tot'.format(combo)]*mask)
        x = np.array(df['{}_unc_tot'.format(combo)]*mask)
    
        x = x[np.isfinite(x)]
        p5.append(np.round(np.percentile(x, 5),2))
        p95.append(np.round(np.percentile(x, 95),2))
    
    print('unc:')
    print('95%: {} to {}\n'.format(np.min(p5),np.max(p95)))
    
    if t0>2002:
        combo = 'JPL'
        print(stats(combo))
    combo ='IMB+WGP'
    print(stats(combo))
#%% unc percentages
periods  = [(2003,2016),
              # (1993,2017)
            ]
for period in periods:
    t0,t1 = period
    print('\nFrom {} to {}'.format(t0,t1))
    path='/Volumes/LaCie_NIOZ/data/barystatic/revisions/results_final/{}-{}/'.format(t0,t1)

    file = 'OM_reconstructions_{}-{}.p'.format(t0,t1)

    df = pd.read_pickle(path+file)

    if t0>2002:
        combos = ['JPL','CSR',
                      'IMB+WGP',
                       'IMB+GWB+ZMP',
                       'UCI+WGP','UCI+GWB+ZMP',
                      ]
    else:
        combos = [#'JPL','CSR',
                      'IMB+WGP',
                       'IMB+GWB+ZMP',
                       'UCI+WGP','UCI+GWB+ZMP',
                      ]   
    unc_type = ['temporal','spatial','intrinsic']
    percs = np.zeros((len(combos),len(unc_type)))
    for i, combo in enumerate(combos):
        tot = np.nanmean(np.array(df['{}_unc_tot'.format(combo)]*mask))
        uncs = np.zeros((len(unc_type)))
        for j,u in enumerate(unc_type):
            uncs[j] = np.nanmean(np.array(df['{}_unc_{}'.format(combo,u)]*mask))
        unc = np.sum(uncs)
        intr = np.nanmean(np.array(df['{}_unc_intrinsic'.format(combo)]*mask))
        spat = np.nanmean(np.array(df['{}_unc_spatial'.format(combo)]*mask))
        temp = np.nanmean(np.array(df['{}_unc_temporal'.format(combo)]*mask))
        unc = np.array(intr+spat+temp)
        percs[i] = np.array((uncs/unc) * 100)
        for j,u in enumerate(unc_type):
            
            percs[i,j] = np.array(( np.nanmean(np.array(df['{}_unc_{}'.format(combo,u)]*mask)) * 100)/unc)
        # percs[i,0] = np.array((temp * 100) /unc)
        # percs[i,1] = np.array((spat * 100) /unc)
        # percs[i,2] = np.array((intr * 100) /unc)
        
    percs[percs==0]=np.nan
    avg_perc = np.nanmean(percs,axis=0)
    print('average uncertainty percentage is:')
    for j,u in enumerate(unc_type):
        print('{}: {:.0f}%'.format(u,avg_perc[j]))
#%% coastal examples
import xarray as xr
#% %
def from_360_to_180(lon_in):
    # given a longitude array that goes from 0 to 360, 
    # returns an array that goes from -180 to 180
    lon_out=np.copy(lon_in)
    for i,ilon in enumerate(lon_in):
        if ilon > 180: 
            lon_out[i]=ilon-360
    return lon_out
#%%
path='/Volumes/LaCie_NIOZ/data/barystatic/revisions/results_final/{}-{}/'.format(t0,t1)

file = 'coastal_examples_10_{}-{}.p'.format(t0,t1)

df = pd.read_pickle(path+file)
t0=2003;t1=2016

t0=2003;t1=2016
path = '/Volumes/LaCie_NIOZ/data/barystatic/revisions/results_final/{}-{}/'.format(t0,t1)
# path = '/Users/ccamargo/Desktop/'
ds=xr.open_dataset(path+'final_dataset_{}-{}.nc'.format(t0,t1))
t0_a=1993
t1_a=2017
path_a = '/Volumes/LaCie_NIOZ/data/barystatic/revisions/results_final/{}-{}/'.format(t0_a,t1_a)
# path = '/Users/ccamargo/Desktop/'
ds_a=xr.open_dataset(path_a+'final_dataset_{}-{}.nc'.format(t0_a,t1_a))



df= pd.read_pickle(path_a+'coastal_examples_10_{}-{}.p'.format(t0_a,t1_a))
df['label'] = ['Cape Town \n(ZA)','Jakarta \n(ID)', 'Lima \n(PE)', 'Mogadishu \n(SO)', 'Rio de Janeiro \n(BR)',
              'Rotterdam \n(NL)', 'Sydney \n(AU)', 'Tokyo \n(JP)', 'Vancouver \n(CA)', 'Washington \n(US)']
df['lon2']=from_360_to_180(df['lon'])
# df=df.sort_values(by=['lat','lon2'],ascending=[False,True]).reset_index()
df= df.sort_values(by=['lon2','lat'],ascending=[True,False]).reset_index()
df=df.drop('index',axis=1)

df_global=pd.read_pickle(path_a+'OM_reconstructions_{}-{}.p'.format(t0_a,t1_a))

name='IMB+WGP_unc'


df_sel=pd.DataFrame({'{}_AIS'.format(name):np.array(df['{}_AIS'.format(name)]),
                     '{}_GIS'.format(name):np.array(df['{}_GIS'.format(name)]),
                     '{}_GLA'.format(name):np.array(df['{}_GLA'.format(name)]),
                     '{}_LWS'.format(name):np.array(df['{}_LWS'.format(name)]),
                     
                     }
                    )
total_unc = np.array(df['{}_tot'.format(name)])
name_trend = name.split('_')[0]+'_trend_tot'

lon=np.array(df_global['lon']).reshape(180,360)[0,:]
lon[-1]=360
lat=np.array(df_global['lat']).reshape(180,360)[:,0]

# path='/Volumes/LaCie_NIOZ/data/barystatic/results/{}-{}/'.format(t0,t1)
# ds=xr.open_dataset(path+'final_dataset_{}-{}.nc'.format(t0,t1))
# ds=ds.sortby('lat',ascending=True)
dataset=name.split('_')[0]
dataset0=dataset.split('+')[0]
dataset1=dataset.split('+')[1]
names=[name for name in np.array(ds.name) if name.split('_')[1]==dataset0 or 
                                             name.split('_')[1]==dataset1]

da=ds.sel(name=names)
locations = df[['lon','lat']]
llon,llat=np.meshgrid(np.array(ds.lon),np.array(-ds.lat))
# AIS
i=0
df_AIS = pd.DataFrame({ '{}'.format(unc):np.array(da.uncs.sel(name=names[i])[ind]).flatten() for ind,unc in enumerate(np.array(ds.unc_type))})
df_AIS['lon']=llon.flatten()
df_AIS['lat']=llat.flatten()
df_AIS=pd.merge(locations,df_AIS,on=['lon','lat'])
df_AIS=df_AIS.drop(['lat','lon'],axis=1)

i=i+1
df_GIS = pd.DataFrame({'{}'.format(unc):np.array(da.uncs.sel(name=names[i])[ind]).flatten() for ind,unc in enumerate(np.array(ds.unc_type))})
df_GIS['lon']=llon.flatten()
df_GIS['lat']=llat.flatten()
df_GIS=pd.merge(locations,df_GIS,on=['lon','lat'])
df_GIS=df_GIS.drop(['lat','lon'],axis=1)
i=i+1
df_GLA = pd.DataFrame({'{}'.format(unc):np.array(da.uncs.sel(name=names[i])[ind]).flatten() for ind,unc in enumerate(np.array(ds.unc_type))})
df_GLA['lon']=llon.flatten()
df_GLA['lat']=llat.flatten()
df_GLA=pd.merge(locations,df_GLA,on=['lon','lat'])
df_GLA=df_GLA.drop(['lat','lon'],axis=1)
i=i+1
df_LWS = pd.DataFrame({'{}'.format(unc):np.array(da.uncs.sel(name=names[i])[ind]).flatten() for ind,unc in enumerate(np.array(ds.unc_type))})
df_LWS['lon']=llon.flatten()
df_LWS['lat']=llat.flatten()
df_LWS=pd.merge(locations,df_LWS,on=['lon','lat'])
df_LWS=df_LWS.drop(['lat','lon'],axis=1)
#%%
df_all = round(df_LWS+df_AIS+df_GIS+df_GLA,2)
total = df_all.sum(axis=1)
al = np.array(df_all)
tot = np.array(total)
np.round(al[:,0]*100/tot)
np.round(al[:,1]*100/tot)
np.round(al[:,2]*100/tot)
#%%
regs=['LWS','AIS','GIS','GLA']
i=0
for df in [df_LWS,df_AIS,df_GIS,df_GLA]:
    print(regs[i],'')
    tg = np.array(df)
    print('temporal:\n')
    print(np.round(tg[:,0]*100/tot))
    print('spatial:\n')
    print(np.round(tg[:,1]*100/tot))
    print('intrinsic:\n')
    print(np.round(tg[:,2]*100/tot))
    print('')
    i=i+1
    
    

#%%
mask=xr.open_dataset('/Volumes/LaCie_NIOZ/PhD/Barystatic/GRACE/mascons/JPL/global/LAND_MASK.CRI_360x180.nc')
# mask.land_mask.plot()
mask=mask.sortby('lat',ascending=False)
oceanmask=np.array(mask.land_mask)
oceanmask[oceanmask>0]=np.nan
oceanmask[oceanmask==0]=1
# plt.pcolor(oceanmask)
da=ds.uncs[2]
names =['AIS_IMB','AIS_JPL',
        'GIS_IMB','GIS_JPL',
        'GLA_JPL',
        'LWS_JPL']
#%%
for name in names:
    da2 = da.sel(name=name)
    data = np.array(da2.data * oceanmask)
    print('{}: \nmin: {:.2f}\nmax: {:.2f}\n mean:{:.2f}'.format(name,np.nanmin(data),
                                                    np.nanmax(data),np.nanmean(data)))
