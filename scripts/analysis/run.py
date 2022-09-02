#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 11:40:29 2022

@author: ccamargo
"""

import os
def run(script_list):
    for script in script_list:
        print('Running {}'.format(script))
        os.system('python {}/{}'.format(path,script))
        
def run_nb(script_list):
    for script in script_list:
        print('Running {}'.format(script))
        os.system('jupyter nbconvert --to notebook --execute {}/{}'.format(path,script))

path = '/Users/ccamargo/Desktop/revisions/m21/'
path = '/Users/ccamargo/Documents/OM_scripts/revisions/m21/'
flist = sorted([file for file in os.listdir(path) if file.endswith('.py')])
if 'run.py' in flist:
    flist.remove('run.py')
# flist = [file for file in flist if int(file.split('-')[0][0])!=4]
# flist = [file for file in flist if int(file.split('-')[0][0])!=3]

run(flist)

# flist = sorted([file for file in os.listdir(path) if file.endswith('.ipynb')])
run_nb(path+'OM_figures-reviews-v4.ipynb')
