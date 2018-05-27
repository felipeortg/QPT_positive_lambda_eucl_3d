#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2018-04-26 17:09:24
# @Author  : Felipe Ortega (felipeortegagama@gmail.com)
# @Link    : ${link}
# @Version : $Id$

import os
import sys
import numpy as np
import csv
from config import *

# Modify this for different data dates
# Begin with _ for better syntax of csv file


def getdata():
    # --------------
    # Get the Eigenvalues

    csv_folder = 'results/csvfiles/'

    # -----------
    # -----------
    # Get data to plot

    print 'Getting processed data from files'

    if onlyintegers:
        prefolder = 'results/pre-plots/Np_'+ str(punc) + '_i_'+ process_date +'/'
    else:
        prefolder = 'results/pre-plots/Np_'+ str(punc) + '_' + process_date +'/'

    name0 = 'Torus'+str(punc)+str(lattice_type)+'_'+str(Kspec)+'_'+str(Ktilde)+'_alps'
    name1 = 'Torus'+str(punc)+str(lattice_type)+'_'+str(Kspec)+'_'+str(Ktilde)
    name2 = 'Torus_de1_'+str(punc)+str(lattice_type)+'_'+str(Kspec)+'_'+str(Ktilde)
    name3 = 'Torus_de2_'+str(punc)+str(lattice_type)+'_'+str(Kspec)+'_'+str(Ktilde)
    name4 = 'Torus_gap_'+str(punc)+str(lattice_type)+'_'+str(Kspec)+'_'+str(Ktilde)
    name5 = 'Torus_order_'+str(punc)+str(lattice_type)+'_'+str(Kspec)+'_'+str(Ktilde)
    name6 = 'Torus_binder_'+str(punc)+str(lattice_type)+'_'+str(Kspec)+'_'+str(Ktilde)
    name7 = 'Torus_susc_'+str(punc)+str(lattice_type)+'_'+str(Kspec)+'_'+str(Ktilde)


    alps = np.load(prefolder + name0 +'.npy')
    Groundstate = np.load(prefolder + name1 +'.npy')
    d1zero = np.load(prefolder + name2 +'.npy')
    d2zero, alpcrit = np.load(prefolder + name3 +'.npy')
    energygap = np.load(prefolder + name4 +'.npy')
    order = np.load(prefolder + name5 +'.npy')
    binder = np.load(prefolder + name6 +'.npy')
    susc = np.load(prefolder + name7 +'.npy')

    print 'Saving CSVs'


    # Save the gap and ground state values
    if not os.path.exists(csv_folder):
        os.makedirs(csv_folder)


    if onlyintegers:
        filename = csv_folder + 'Np_'+ str(int(punc)) +'_i_'+ process_date + '.csv'
    else:
        filename = csv_folder + 'Np_'+ str(int(punc)) +'_'+ process_date + '.csv'
    
    with open(filename, 'w') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['param','Ground Energy', 'Gap energy', 'Sec Der of E0/punc', 'Order param', 'Binder', 'Susceptibility'])

        for nn, ens in enumerate(energygap[1]):
            writer.writerow([round(alps[nn],4),Groundstate[nn][0],ens,d2zero[nn],order[nn],binder[nn],susc[nn]])



if __name__ == '__main__':

    if len(sys.argv) > 1:
        raise ValueError('Invalid number of arguments')
    else:
        getdata()
