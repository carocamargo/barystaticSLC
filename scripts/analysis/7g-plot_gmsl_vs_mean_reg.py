#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 23 10:03:26 2021

@author: ccamargo
"""


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 22 15:14:59 2021

@author: ccamargo
"""



import pandas as pd
import numpy as np 
import xarray as xr
import sys
sys.path.append("/Users/ccamargo/Documents/py_scripts/")
import utils_SLE_v2 as sle


def quadrsum(X):
    if len(X.shape)==3:
        Y = np.zeros((X.shape[1],X.shape[2]))
    elif len(X.shape)==2:
        Y = np.zeros((X.shape[1]))
    else:
        Y = 0
    for i in range(X.shape[0]):
        Y = Y + (X[i]*X[i])
    Y = np.sqrt(Y)
    return Y
import matplotlib.pyplot as plt

col_dict={1:"black", # WN
          2:"crimson", # PL
          3:"skyblue", # PLWN
          4:"orange", # AR1
          5:"teal", # Ar5
          6:"darkmagenta", # AR9
          
          7:"skyblue", # ARf
          8:"crimson" # GGM
          }



#%% mask
mask=xr.open_dataset('/Users/ccamargo/Documents/PhD/Barystatic/GRACE/mascons/JPL/global/LAND_MASK.CRI_360x180.nc')
# mask.land_mask.plot()
mask=mask.sortby('lat',ascending=False)
oceanmask=np.array(mask.land_mask)
oceanmask[oceanmask>0]=np.nan
oceanmask[oceanmask==0]=1
#%%
t0=2003
t1=2016
# path = '/Volumes/LaCie_NIOZ/data/barystatic/results/{}-{}/'.format(t0,t1)
# # path = '/Users/ccamargo/Desktop/'
# ds=xr.open_dataset(path+'final_dataset_{}-{}.nc'.format(t0,t1))

path = '/Volumes/LaCie_NIOZ/data/barystatic/'
df=pd.read_pickle(path+'global/final_combos-{}-{}_v2.p'.format(t0,t1))

df_reg = pd.read_pickle(path+'results/{}-{}/OM_reconstructions_{}-{}.p'.format(t0,t1,t0,t1))

dimlat=180;dimlon=360
lon=np.array(df_reg['lon']).reshape(dimlat,dimlon)[0,:]
lat=np.array(df_reg['lat']).reshape(dimlat,dimlon)[:,0]

reg_trend=np.zeros((len(df['combinations'])))
# reg_trend_ols=np.full_like(reg_trend,0)
reg_unc=np.full_like(reg_trend,0)

for i,combo in enumerate(np.array(df['combinations'])):
    reg_trend[i],_,_=sle.reg_to_glb(np.array(df_reg[combo+'_trend_tot']).reshape(dimlat,dimlon),lat,lon)
    reg_unc[i],_,_=sle.reg_to_glb(np.array(df_reg[combo+'_unc_tot']).reshape(dimlat,dimlon),lat,lon)
    # reg_trend_ols[i],_,_=sle.reg_to_glb(np.array(df_reg_ols[combo+'_trend_tot']).reshape(dimlat,dimlon),lat,lon)

df['reg_trend_nm'] = reg_trend
# df['reg_trend_ols'] = df_reg_ols['trend']
df['reg_unc'] = reg_unc


#%%
col = ['trend','reg_trend_nm',
       'trend_OLS','trend_OLS_reg',
       'total_unc','reg_unc',
       # 'temp_unc',
       ]
labels = ['gmsl trend (NM)', 'mean(regional trend(NM))', 
          'gmsl trend (OLS)','mean(regional trend(OLS))', 
          'gmsl unc', 'mean (regional unc)' ,
          ]
fig=plt.figure(figsize=(10,5),dpi=100,facecolor='w')

syb=['s','s','d','d',"^",'^','v','*','D','X','X','X']
ax=plt.subplot(111)
for i,c in enumerate(col):
    df.plot.scatter(x='combinations',y=c,ax=ax,
                    marker=syb[i],
                    s=80,
                    color=col_dict[i+1],
                    label=labels[i])
plt.legend(loc='upper left',ncol=int(len(labels)/2),fontsize=12)
plt.ylim([0,3])
plt.xticks(np.arange(0,len(df['combinations'])),
           df['combinations'],
           rotation=75,fontsize=12)
plt.yticks(fontsize=12)
plt.ylabel('mm/year',fontsize=12)
plt.xlabel('combinations',fontsize=12)
plt.grid()
plt.title('Global mean sea-level trend')
plt.show()
plt.close()
#%%
#%%
col = ['trend','reg_trend_nm',
       # 'trend_OLS','reg_trend_ols',
       'total_unc','reg_unc',
       # 'temp_unc',
       ]
labels = ['GMSL trend (NM)', 'mean (regional trend)', 
          # 'gmsl trend (OLS)','mean(regional trend(OLS))', 
          'GMSL uncertainty','mean (regional uncertainty)' ,
           ]
fig=plt.figure(figsize=(10,5),dpi=100,facecolor='w')

syb=['s','s',
     # 'd','d',
     "^",'^','v','*','D','X','X','X']
ax=plt.subplot(111)
for i,c in enumerate(col):
    df.plot.scatter(x='combinations',y=c,ax=ax,
                    marker=syb[i],
                    s=80,
                    color=col_dict[i+1],
                    label=labels[i])
plt.legend(loc='upper left',ncol=int(len(labels)/2),fontsize=12)
plt.ylim([0,2.5])
plt.xticks(np.arange(0,len(df['combinations'])),
           df['combinations'],
           rotation=75,fontsize=12)
plt.yticks(fontsize=12)
plt.ylabel('mm/year',fontsize=12)
plt.xlabel('combinations',fontsize=12)
plt.grid()
plt.show()
#%%
path_to_figures = '/Users/ccamargo/Documents/PhD/Barystatic/manuscript/latex/figures_v3/overleaf/'
fig.savefig(path_to_figures+'gmsl_{}-{}.png'.format(t0,t1),dpi=300,format='png')
plt.close()
#%%
col = ['trend','reg_trend_nm',
       # 'trend_OLS',
       'trend_OLS_reg',
       'total_unc','reg_unc',
       # 'temp_unc',
       ]
labels = ['gmsl trend (NM)', 'mean(regional trend(NM))', 
          # 'gmsl trend (OLS)',
          'mean(regional trend(OLS))', 
          'gmsl unc','mean (regional unc)' ,
           ]
fig=plt.figure(figsize=(10,5),dpi=100,facecolor='w')

syb=['s','s',
     # 'd',
     'd',
     "^",'^','v','*','D','X','X','X']
ax=plt.subplot(111)
for i,c in enumerate(col):
    df.plot.scatter(x='combinations',y=c,ax=ax,
                    marker=syb[i],
                    s=80,
                    color=col_dict[i+1],
                    label=labels[i])
plt.legend(loc='upper left',ncol=int(len(labels)/2),fontsize=12)
plt.ylim([0,3.5])
plt.xticks(np.arange(0,len(df['combinations'])),
           df['combinations'],
           rotation=75,fontsize=12)
plt.yticks(fontsize=12)
plt.ylabel('mm/year',fontsize=12)
plt.xlabel('combinations',fontsize=12)
plt.grid()
plt.show()
plt.close()
