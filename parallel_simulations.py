

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 26 13:56:56 2022

@author: amelie aussel
"""


print('Starting')

import os


import matplotlib
import matplotlib.pyplot as plt
import time
#matplotlib.use('agg')
plt.switch_backend('agg')
import sys

from brian2 import *

from joblib import Parallel, delayed
import multiprocessing
import os

from model_files.FEF_and_LIP_single_simulation import *

os.environ['MKL_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['MKL_DYNAMIC'] = 'FALSE'

import time
import ntpath
from itertools import *




path=""
if os.name == 'nt':
    path=os.path.join(ntpath.dirname(os.path.abspath(__file__)),"results_"+str(datetime.datetime.now()).replace(':','-'))
else :
    path="/results_"+str(datetime.datetime.now())

os.mkdir(path)

#list of target presentation times :
liste_target_time=[300*msecond,400*msecond,500*msecond,600*msecond,700*msecond,800*msecond,900*msecond,1000*msecond,1100*msecond,1200*msecond,1300*msecond,1400*msecond,1500*msecond,1600*msecond,1700*msecond]
liste_target_time+=[350*msecond,450*msecond,550*msecond,650*msecond,750*msecond,850*msecond,950*msecond,1050*msecond,1150*msecond,1250*msecond,1350*msecond,1450*msecond,1550*msecond,1650*msecond]
liste_target_time+=[325*msecond,425*msecond,525*msecond,625*msecond,725*msecond,825*msecond,925*msecond,1025*msecond,1125*msecond,1225*msecond,1325*msecond,1425*msecond,1525*msecond,1625*msecond,1725*msecond]
liste_target_time+=[375*msecond,475*msecond,575*msecond,675*msecond,775*msecond,875*msecond,975*msecond,1075*msecond,1175*msecond,1275*msecond,1375*msecond,1475*msecond,1575*msecond,1675*msecond]
# N simulations will be ran for each possible simulation time :
N=50


#When varying tFS and tSOM in the paper, less simulations were performed
#N=20
#liste_target_time=[300*msecond,400*msecond,500*msecond,600*msecond,700*msecond,800*msecond,900*msecond,1000*msecond,1100*msecond,1200*msecond,1300*msecond,1400*msecond,1500*msecond,1600*msecond,1700*msecond]
#liste_target_time+=[350*msecond,450*msecond,550*msecond,650*msecond,750*msecond,850*msecond,950*msecond,1050*msecond,1150*msecond,1250*msecond,1350*msecond,1450*msecond,1550*msecond,1650*msecond]


liste_simus=[]
for t in liste_target_time:
    liste_simus+=[t]*N

#Setting up the lists of SOM and FS inhibition decay times to use :
liste_t_SOM=[20*ms]
liste_t_FS=[5*ms]

#Other parameters (fixed across all simulations):
theta_phase='mixed' #theta phases to simulate (good, bad or mixed)
gLIP_FEFv=0.015*msiemens * cm **-2 #LIP->FEF visual module synapse conductance :
target_presentation='True' #'True' if target is presented, 'False' otherwise
runtime=2*second #simulation duration


liste_simus=[[i,j,k] for i in liste_simus for j in liste_t_SOM for k in liste_t_FS]
liste_simus=[[liste_simus[i][0],i,liste_simus[i][1],liste_simus[i][2],theta_phase,gLIP_FEFv,target_presentation,runtime] for i in range(len(liste_simus))]


#setting the number of cores to used (all cpus by default)
num_cores = multiprocessing.cpu_count()
if os.name == 'nt':
    num_cores=-3 #using all cpus on a windows does not work for an unknown reason

Parallel(n_jobs=num_cores)(delayed(FEF_and_LIP)(simu,path) for simu in liste_simus)


clear_cache('cython')
