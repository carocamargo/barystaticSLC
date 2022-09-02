#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 18 16:34:47 2022

@author: ccamargo
"""
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import seaborn as sns
import sys
sys.path.append("/Users/ccamargo/Documents/py_scripts/")
import utils_SL as sl 
from matplotlib.gridspec import GridSpec
#%%
periods =[#(2005,2016),(1993,2018),
          (2003,2017)
          ]
for period in periods:
    t0=period[0]
    t1=period[1]-1
    
    pwd = '/Volumes/LaCie_NIOZ/data/barystatic/revisions/hector/'
    
    pwd = pwd+'{}-{}/'.format(t0,t1)

    flist=sl.get_filelist(pwd+'ranking/regional/','*.nc')
    name=[None]*len(flist)
    ifile=11;file=flist[ifile]

    dimlat=180;dimlon=360
    nrows=4
    ncols=4
    wi=20;hi=15
    for score in [#'aic','bic',
                  'bic_tp'
                  ]:
        fig = plt.figure(figsize=(wi, hi), 
                     facecolor='w')
    
        gs = GridSpec(nrows, ncols)
        for ifile, file in enumerate(flist):
            ds=xr.open_dataset(file)
            name = file.split('/')[-1].split('.nc')[0]
            data = {}
            for i, nm in enumerate(np.array(ds.nm)):
                datab = np.array(ds[score][i])
                datab = datab[np.isfinite(datab)]
                datab = datab[datab>0]
                data[nm] = datab
            
            ax1 = plt.subplot(gs[ifile])
            # if name in ['AIS_IMB','AIS_UCI','GIS_IMB','GIS_UCI']:
            #     mode='layer'
            # else:
            #     mode='stack'
            sns.histplot(data,
                          multiple='stack',
                         # stat='density'
                         )
            plt.title(name)
        plt.show()