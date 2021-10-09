#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 22 15:14:59 2021

@author: ccamargo
"""



import pandas as pd
import numpy as np 
import xarray as xr
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
#%%
# t0=2005
# t1=2015

periods =[ 
            # (2005,2016),
            # (1993,2018),
            (1993,2017),
            (2003,2017)
          ]


combos = [['JPL'],
            ['CSR'],
           ['IMB','WGP'],
          ['IMB','GWB','ZMP'],['UCI','WGP'],['UCI','GWB','ZMP']
         ]
reconstr =['JPL',
            'CSR',
            'IMB+WGP',
           'IMB+GWB+ZMP','UCI+WGP','UCI+GWB+ZMP'
          ]

period=periods[1]
for period in periods:
    t0=period[0]
    t1=period[1]-1
    print(period)
    path = '/Volumes/LaCie_NIOZ/data/barystatic/global/'
    #% 
    ds_temporal = xr.open_dataset(path+'hector/ALL-ocean_mean_8_NM-{}-{}_v2.nc'.format(t0,t1))
    ds_temporal = xr.open_dataset(path+'hector/ALL-ocean_mean_9_NM-{}-{}_v2.nc'.format(t0,t1)) # OLS
    ic_idx=-1    
    df=pd.DataFrame(np.array(ds_temporal.name),columns=['dataset'])
    df['trend']=np.array(ds_temporal.best_tr[:,ic_idx])
    df['temporal_unc']=np.array(ds_temporal.best_unc[:,ic_idx])
    df['trend_OLS']=np.array(ds_temporal.trend[:,-1])
    df_reg_ols = pd.read_pickle(path+'OLS_SLF_reg_mean-{}-{}.p'.format(t0,t1))
    df_reg_ols = df_reg_ols.rename(columns={'trend':'trend_OLS_reg'})
    df = df.merge(df_reg_ols,on='dataset')
    
    df_intrinsic = pd.read_pickle(path+'intrinsic_OLS_trend_{}-{}_v2.p'.format(t0,t1))
    df_intrinsic=df_intrinsic.rename(columns={'trend':'intrinsic_unc'})
    df = df.merge(df_intrinsic,on='dataset',how='left')
    df_structural = pd.read_pickle(path+'structural_HEC_trend_{}-{}_v2.p'.format(t0,t1))   
    df['structural_unc'] = [df_structural.loc[reg.split('_')[0]][0] for reg in np.array(df['dataset'])]
    df=df.fillna(0)
    df['region'] = [reg.split('_')[0] for reg in np.array(df['dataset'])]
    X= np.array(df[['temporal_unc','intrinsic_unc','structural_unc']]).T
    df['sum_unc']= quadrsum(X)
    X= np.array(df[['temporal_unc','intrinsic_unc']]).T
    df['sum_unc_temp_intr']= quadrsum(X)
    
    df.to_pickle(path+'final_dataset-{}-{}_v2.p'.format(t0,t1))
    
    #% %
    df_combo =pd.DataFrame(reconstr,columns=['combinations'])
    df_combo['dataset'] = combos
    # names = [None]*len(combos)
    names=[]
    total_unc = np.zeros((len(combos)))
    temp_unc = np.zeros((len(combos)))
    temp_intr_unc = np.zeros((len(combos)))
    trend=np.zeros((len(combos)))
    trend_ols=np.zeros((len(combos)))
    trend_ols_reg=np.zeros((len(combos)))
    
    for i,combo in enumerate(combos):
        # print(title)
        #% %
        names.append([name for name in np.array(df['dataset']) for comb in combo if name.split('_')[1]==comb])
        df_sel = df.loc[df['dataset'].isin(names[i])].reset_index()
        # X= np.array(df_sel[['temporal_unc','intrinsic_unc','structural_unc']]).T
        # df_sel['sum_by_type']= quadrsum(X)
        total_unc[i] = quadrsum(np.array(df_sel['sum_unc']).T)
        temp_unc[i] = quadrsum(np.array(df_sel['temporal_unc']).T)
        temp_intr_unc[i] = quadrsum(np.array(df_sel['sum_unc_temp_intr']).T)
        trend[i] = np.sum(np.array(df_sel['trend']))
        trend_ols[i] = np.sum(np.array(df_sel['trend_OLS']))
        trend_ols_reg[i] = np.sum(np.array(df_sel['trend_OLS_reg']))
    
    df_combo['trend'] = trend
    df_combo['total_unc']=total_unc
    df_combo['temp_unc']=temp_unc
    df_combo['temp_intr_unc']=temp_intr_unc
    df_combo['trend_OLS'] = trend_ols
    df_combo['trend_OLS_reg'] = trend_ols_reg

    df_combo['names']=names
    for i,combo in enumerate(np.array(df_combo['combinations'])):
        # print(combo)
        print('{}: {:.2f} ± {:.2f}'.format(combo, df_combo.iloc[i]['trend'], 
                                           df_combo.iloc[i]['total_unc']))
    # 
    df_tmp = df[['dataset','trend','sum_unc','region']]
    df_tmp_c= df_combo[['combinations','trend','total_unc','names']]
    df_tmp = df_tmp.rename(columns={'sum_unc':'total_unc'})
    df_tmp_c.rename(columns={'combinations':'dataset','names':'datasets'},inplace=True)
    df2= df_tmp.append(df_tmp_c)
    df2.to_excel(path+'global_mean_table-{}-{}.xlsx'.format(t0,t1))
    #% %
    df_combo.to_pickle(path+'final_combos-{}-{}_v2.p'.format(t0,t1))
