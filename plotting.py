#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2018-05-11 23:02:06
# @Author  : Felipe Ortega (felipeortegagama@gmail.com)
# @Title   : Plot


import matplotlib.pyplot as plt

from config import *


# -----------------
# Matplotlib param
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{physics}')
plt.rc('font', size=14)
plt.rc('font', family='serif')
plt.rc('axes.formatter', useoffset = False)


def plots(save_show=0):

    # --------------
    # Where the plots will be saved
    if onlyintegers:
        save_folder = 'results/plots/Np_'+ str(punc) + '_i_/'
    else:
        save_folder = 'results/plots/Np_'+ str(punc)  +'/'


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

    alps = np.load(prefolder + name0 +'.npy')
    Groundstate = np.load(prefolder + name1 +'.npy')
    d1zero = np.load(prefolder + name2 +'.npy')
    d2zero, alpcrit = np.load(prefolder + name3 +'.npy')
    energygap = np.load(prefolder + name4 +'.npy')
    order = np.load(prefolder + name5 +'.npy')

    TotalEigs = len(Groundstate[0])
    # How many intervals (number of points -1)
    maxin = len(Groundstate) - 1

    print 'Plotting'

    # --------------
    # --------------
    # Plot Spectrum
    fig1 = plt.figure(1,figsize=(7, 5))

    for nn in xrange(TotalEigs):
        plt.plot(alps, [Groundstate[mm][nn] for mm in xrange(maxin+1)],'-')


    # tidy up the figure
    plt.grid(True)
    #plt.legend(loc='upper right')
    plt.title(r'$H_{'+str(punc)+str(lattice_type)+r'} = -\lambda \sum_p B_p(k = '+ str(Kspec)+r') - (1-\lambda) \sum_p \tilde{B}_p(\tilde{k} = '+ str(Ktilde)+r')$')
    plt.xlabel(r'$\lambda$')
    plt.ylabel(r'$\textrm{Energy/puncture}$')
    # plt.axis([0,1,-3.2,0.2])
    # plt.axes().set_aspect(0.05)

    figname1 = 'Torus'+str(punc)+str(lattice_type)+'_'+str(Kspec)+'_'+str(Ktilde)+'.pdf'


    # Title for the rest of the figures
    title_name = r'$H_{'+str(punc)+str(lattice_type)+r'},\, (k = '+ str(Kspec)+r',\tilde{k} = '+ str(Ktilde)+r'),\,\lambda_c \approx '+ str(alpcrit) +r' $'
    # Plot the first derivative

    fig2 = plt.figure(2)
 
    plt.plot(alps[2:-2],d1zero[2:-2])

    plt.grid(True)
    plt.title(title_name)
    plt.ylabel(r'$\frac{\partial E_0}{\partial\lambda}$')
    plt.xlabel(r'$\lambda$')
    fig2.subplots_adjust(left=0.17)

    figname2 = 'Torus_de1_'+str(punc)+str(lattice_type)+'_'+str(Kspec)+'_'+str(Ktilde)+'.pdf'


    # Plot the second derivative

    fig3 = plt.figure(3)

    plt.plot(alps[1:-1],d2zero[1:-1])

    plt.grid(True)
    plt.title(title_name)
    plt.ylabel(r'$\frac{\partial^2E_0}{\partial\lambda^2}$')
    plt.xlabel(r'$\lambda$')
    fig3.subplots_adjust(left=0.17)


    figname3 = 'Torus_de2_'+str(punc)+str(lattice_type)+'_'+str(Kspec)+'_'+str(Ktilde)+'.pdf'

    # Plot the gap

    fig4 = plt.figure(4)

    for nn in xrange(TotalEigs-1):
        plt.plot(alps,energygap[nn+1],'-')

    plt.grid(True)
    plt.title(title_name)
    plt.xlabel(r'$\lambda$')
    plt.ylabel(r'Excitation energies')

    figname4 = 'Torus_gap_'+str(punc)+str(lattice_type)+'_'+str(Kspec)+'_'+str(Ktilde)+'.pdf'

    # Plot the order parameter

    fig5 = plt.figure(5)

    plt.plot(alps,order,'-')

    plt.grid(True)
    plt.title(title_name)
    plt.xlabel(r'$\lambda$')
    if onlyintegers:
        plt.ylabel(r'$\expval{W_{1}}$')
    else:
        plt.ylabel(r'$\expval{W_{1/2}}$')

    figname5 = 'Torus_order_'+str(punc)+str(lattice_type)+'_'+str(Kspec)+'_'+str(Ktilde)+'.pdf'

    if save_show==0:

        if not os.path.exists(save_folder):
            os.makedirs(save_folder)

        fig1.savefig(save_folder + figname1)
        fig2.savefig(save_folder + figname2)
        fig3.savefig(save_folder + figname3)
        fig4.savefig(save_folder + figname4)
        fig5.savefig(save_folder + figname5)
        
    else:

        plt.show()


if __name__ == '__main__':

    if len(sys.argv) > 2:
        raise ValueError('Invalid number of arguments')
    else:
        plots(save_show)







