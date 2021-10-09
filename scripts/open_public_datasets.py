#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  9 09:25:23 2021

Datasets used for figures 

@author: ccamargo
"""

import pandas as pd
import xarray as xr

#% % load data

# periods
t0=2003;t1=2016
t0_a=1993;t1_a=2017

# paths
path = '/Volumes/LaCie_NIOZ/data/barystatic/results/{}-{}/'.format(t0,t1)
path_a = '/Volumes/LaCie_NIOZ/data/barystatic/results/{}-{}/'.format(t0_a,t1_a)

# results
ds=xr.open_dataset(path+'final_dataset_OLS-prop_{}-{}.nc'.format(t0,t1))
ds_a=xr.open_dataset(path_a+'final_dataset_OLS-prop_{}-{}.nc'.format(t0_a,t1_a))

# global reconstruction
df_global=pd.read_pickle(path_a+'OM_reconstructions_OLS-prop_{}-{}.p'.format(t0_a,t1_a))
# coastal examples
df= pd.read_pickle(path_a+'coastal_examples_10-prop_{}-{}.p'.format(t0_a,t1_a))

#% % land mask
mask=xr.open_dataset('/Users/ccamargo/Documents/PhD/Barystatic/GRACE/mascons/JPL/global/LAND_MASK.CRI_360x180.nc')
