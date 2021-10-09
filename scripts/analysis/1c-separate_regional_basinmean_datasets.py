#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  9 09:51:36 2021

@author: ccamargo
"""


import xarray as xr
import numpy as np


# import warnings
import sys
sys.path.append("/Users/ccamargo/Documents/py_scripts/")
import utils_SL as sl 

# import utils_SLE as sle 
import pandas as pd

#%% REGIONAL DATASETS
#######################

path='/Volumes/LaCie_NIOZ/data/barystatic/use/comb/'
da=xr.open_dataset(path+'ALL_datasets_1993-2020_180x360_v3_update.nc')
print(da.name)
# name = [name for name in np.array(da.name) if name.endswith('CSR') or 
#         name.endswith('JPL') or name.endswith('GWB') or name.endswith('gl')]

da=da.sel(name=[ # 'AIS_IMB', 'AIS_R19_basins',
                # 'GLWS_ZMP', 
                'GLWS_WGP_gl', 
                'AIS_300_CSR', 'GIS_300_CSR', 'LWS_CSR', 'GLWS_CSR', 'TCWS_CSR', 
                # 'AIS_CSR', 'GIS_CSR',
                'AIS_300_JPL', 'GIS_300_JPL', 'LWS_JPL', 'GLWS_JPL', 'TCWS_JPL',
                # 'AIS_JPL', 'GIS_JPL', 
                # 'AIS_proj_CSR', 'GIS_proj_CSR', 'AIS_proj_JPL','GIS_proj_JPL', 
                # 'GIS_IMB', 'GIS_M19', 
                'LWS_GWB', 'LWS_WGP_gl',
                
               # 'TCWS_WaterGAP'
               ])
print(da)
# da.attrs['units']='mm of SLC'
da.attrs['script']='land_mass_variations_regional.py'
da.attrs['metadata']='Regional SL variations for GRACE and Hydrological Models'
da.to_netcdf('/Volumes/LaCie_NIOZ/data/barystatic/source_var/regional_v3.nc')


##################
#%% DATASETS WITH BASIN/REGIONAL MEANS

#%% Greenland
dataset='GIS_IMB'
# # IMBIE 2019 Greendland Dataset
# #This spreadsheet contains the IMBIE-2019 datasets for Greenland, 
# #which includes data on the annual rate of change and cumulative change in Greenland’s ice sheet mass, 
# #its surface mass balance and ice discharge anomalies, and their estimated uncertainty. 

# # Sheet 2: equivalent mean global sea level rise (mm/yr) and sea-level rise (mm)
# # and in units of equivalent mean global sea level rise (millimetres per year – sheet 2, columns B, C, F, G, J and K, 
# #and millimetres – sheet 2, columns D, E, H, I, L and M).
df=pd.read_excel('/Volumes/LaCie_NIOZ/data/barystatic/original/GRE_IMBIE.xlsx',
                  #sheet_name=[0,1] # load first and second sheet as a dict of df
                  sheet_name=1 # sheet 2
                  )

df = df.set_index('Year')
df = df.dropna() # From 1992 onwards

gre_regions = ['GIS','SMB','DYN']
cols=[2,6,10]
slc=np.zeros((len(gre_regions),len(df)))
unc=np.full_like(slc,0)

for i,col in enumerate(cols):
    print(df.columns[col])
    slc[i,:]=np.array(df[df.columns[col]])
    unc[i,:]=np.array(df[df.columns[col+1]])

da=xr.Dataset(data_vars={'{}_slc'.format(dataset):(('reg_{}'.format(dataset),'time'),slc),
                        '{}_unc'.format(dataset):(('reg_{}'.format(dataset),'time'),unc),                           
                        },
                          coords={
                                  'time':np.array(df.index),
                                  'reg_{}'.format(dataset):gre_regions
                                  })

time_dim=len(df)
#%% get ANT dataset
# IMBIE 2018 Antartic Dataset
#This spreadsheet contains the IMBIE-2018 datasets for Antarctica (worksheet 1),
#Antarctic Peninsula (worksheet 2), East Antarctica (worksheet 3) and West Antarctica (worksheet 4). 
#Each worksheet includes data on monthly cumulative ice sheet mass changes and their estimated uncertainty. 
#The data are expressed in units of mass (Gigatons – columns B and C) 
#and in units of equivalent mean global sea level rise (millimetres – columns D and E).

dataset='AIS_IMB'
file='/Volumes/LaCie_NIOZ/data/barystatic/original/ANT_IMBIE.xlsx'
regions = ['AIS','AP','EAIS','WAIS']
local_dim = len(pd.read_excel(file))
unc=np.zeros((len(regions),time_dim))
slc = np.full_like(unc,0)

for sheet, reg in enumerate(regions):
    df=pd.read_excel(file,
                  #sheet_name=[0,1] # load first and second sheet as a dict of df
                  sheet_name=sheet # sheet 1
                  )
    unc[sheet,0:local_dim] = np.array(df['Cumulative sea level contribution uncertainty (mm)'])
    slc[sheet,0:local_dim] = np.array(df['Cumulative sea level contribution (mm)'])

da['{}_slc'.format(dataset)]=(('reg_{}'.format(dataset),'time'),slc)
da['{}_unc'.format(dataset)]=(('reg_{}'.format(dataset),'time'),unc)

da = da.assign_coords({'reg_{}'.format(dataset):regions})

#%% GRE M19 - UCI
file='/Volumes/LaCie_NIOZ/data/barystatic/original/GRE_mouginot2019.xlsx'
dataset='GIS_UCI'
# yearly datasets 
sheet='TOTALloss'
df=pd.read_excel('/Volumes/LaCie_NIOZ/data/barystatic/original/GRE_mouginot2019.xlsx',
                     sheet_name=sheet,# open TOTAL Loss
                     sep='\s+',
                     header=1
                     )
df=df.rename(columns={df.columns[0]:"region"})
df=df.set_index('region')
regions = np.array(df.index)
# remove before 1979: (start year of ANT Rignot)
df = df.drop([year for year in df.columns if year<1979],axis =1 )
years=np.array(df.columns)
slc=np.zeros((len(regions),len(years)))

for ireg,reg in enumerate(regions):
    slc[ireg,:] = np.array(df.loc[reg]/362) # from Gt to mm
    
da['{}_slc'.format(dataset)]=(('reg_{}'.format(dataset),'year'),slc)
da = da.assign_coords({'reg_{}'.format(dataset):regions})
da = da.assign_coords({'year':years})


#%%ANT RIGNOT
file='/Volumes/LaCie_NIOZ/data/barystatic/original/ANT_rignot2019.xlsx'
dataset='AIS_UCI'

# REGIONS:
# df=pd.read_excel(file,
#                      sheet_name=1,
#                      sep='\s+',
#                      header=2
#                      )
# df['Glacier name'][4]='AIS'
# regions = np.array(df['Glacier name'])
# df = df.set_index('Glacier name')
# slc=np.zeros((len(regions),len(years)))
# years = [year for year in df.columns if type(year)==int]

# discharge=df[years]
# smb = df['SMB 1979-2008']
# # df[1979]-df['SMB 1979-2008']
# for ireg, reg in enumerate(regions):
#     mb = np.array(df['SMB 1979-2008'].loc[reg]) - np.array(df[years].loc[reg])
    
#     for iyear in range(1,len(years)):
#         mb[iyear] = mb[iyear] + mb[iyear-1]
#     slc[ireg,0:len(years)]=mb/362 

# da['{}_slc'.format(dataset)]=(('reg_{}'.format(dataset),'year'),slc)
# da = da.assign_coords({'reg_{}'.format(dataset):regions})

#% BASINS
df=pd.read_excel(file,
                     sheet_name=2,
                      sep='\s+',
                     header=2
                     )


df = df.groupby(by=['Basin']).sum()
basins = ["H-H'",'F-G',"E-E'","D-D'","C'-D","B-C","A-A'",
          'J"-K',"G-H","D'-E","A'-B","C-C'","K-A",'J-J"','I"J','I-I"',"H'-I","E'-F"]
slc=np.zeros((len(basins),len(years)))
years = [year for year in df.columns if type(year)==int]

#% % Take the sum per basin:
df_basin = df.groupby(by=['Basin']).sum()

#% % make a data frame with the codes of our mask and what they mean: 
df_codes=pd.DataFrame(basins,columns=['Basin'])

codes = np.arange(1,len(basins)+1)
df_codes['codes']=codes
df_codes=df_codes.set_index('Basin')

# now merge with the dataframe per basin 
df =  pd.concat([df_codes,df_basin],axis=1)

discharge=df[years]
smb = df['SMB 1979-2008']
# df[1979]-df['SMB 1979-2008']
for ireg, reg in enumerate(basins):
    mb = np.array(df['SMB 1979-2008'].loc[reg]) - np.array(df[years].loc[reg])
    
    for iyear in range(1,len(years)):
        mb[iyear] = mb[iyear] + mb[iyear-1]
    slc[ireg,0:len(years)]=mb/362 

da['{}_slc'.format(dataset)]=(('reg_{}'.format(dataset),'year'),slc)
da = da.assign_coords({'reg_{}'.format(dataset):basins})

#%% ZMP
#% %
dataset='GLA_ZMP'
flist=sl.get_filelist('/Volumes/LaCie_NIOZ/data/barystatic/original/GLA_Zemp2019-v1.0/',ext='*.csv')

flist2=[None]*(len(flist)-1)
flist2[0:9]=flist[11:20]
flist2[9:len(flist2)]=flist[1:11]
flist2


#% % 
mchange_per_glacier=np.zeros((len(flist2),len(da.year)))
area_per_glacier=np.zeros((len(flist2),len(da.year)))

for i,f in enumerate(flist2):
    df=pd.read_csv(f, header=26)
    df.head()
    df = df.set_index(df.Year)
    df = df.drop([year for year in df.index if year<1979],axis =0 )
    data=np.array(df[" AW_mwe"]) # mm 
    # data=np.array(df[" AW_Gt"])/362.5 # Gt to mm
    area=np.array(df[" Area_AW_ref_km2"])
    # area=np.array(df[" Area_LW_km2"])

    data=data[0:len(data)]
    area=area[0:len(area)]
    mchange_per_glacier[i,0:len(data)]=data
    area_per_glacier[i,0:len(area)]=area
#% %
gla_regions = ['{:02}_{}'.format(int(f.split('.')[1].split('_')[-2]),f.split('.')[1].split('_')[-1]) 
               for f in flist2]

da['{}_slc'.format(dataset)]=(('reg_{}'.format(dataset),'year'),mchange_per_glacier)
da = da.assign_coords({'reg_{}'.format(dataset):gla_regions})


#%% save da
da.attrs['units']='mm of SLC'
da.attrs['script']='land_mass_variations_means.py'
da.attrs['metadata']='SL variations for AIS, GIS and Glaciers, from IMBIE, UCI (Mouginot and Rignot 2019) and ZEMP 2019'
da.attrs['means']='Values in the regions that they are given (not regionally redistributed yet)'
da.to_netcdf('/Volumes/LaCie_NIOZ/data/barystatic/source_var/means_v2.nc')
