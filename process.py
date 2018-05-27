#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2018-04-19 16:22:25
# @Author  : Felipe Ortega (felipeortegagama@gmail.com)
# Title    : Processing


import scipy.sparse as spar
from config import *

def process():

    # --------------
    # Where the results will be saved
    if onlyintegers:
        save_folder = 'results/pre-plots/Np_'+ str(punc) + '_i_'+ process_date +'/'
    else:
        save_folder = 'results/pre-plots/Np_'+ str(punc) + '_' + process_date +'/'

    # -----------
    # -----------
    # Get order parameter matrix

    print 'Getting B-half ME from files'

    if onlyintegers:
        matrixfolder = 'results/matrices/Np_'+ str(punc) + '_i_'+ B_date +'/'
    else:
        matrixfolder = 'results/matrices/Np_'+ str(punc) + '_' + B_date +'/'

    # Get the B-half
    filename = 'LT_' + str(lattice_type) + '_Trp_' + str(TOPRIGHT) + '_K_'+ str(Kspec)

    mattype = np.load(matrixfolder+filename+'.npy')[0]

    if mattype == 'dense':

        filename +=  '_B-half'

        Bfmathalf = np.load(matrixfolder+filename+'.npy')

        Bfmathalfsq = [ Bfmathalf[ff].dot(Bfmathalf[ff])
            for ff in xrange(punc)]

        # Size of the Hilbert space
        NT = len(Bfmathalf[0])

    elif mattype == 'sparse':

        filename +=  '_B-half'

        # Size of the Hilbert space
        Bfmathalfdata,rowssparse,colssparse,NT = np.load(matrixfolder+filename+'.npy')
        

        Bfmathalf = [
            spar.coo_matrix(
                (Bfmathalfdata[ff],
                    (rowssparse[ff],
                    colssparse[ff])
                ),
            shape = (NT,NT)).tocsr()

            for ff in xrange(punc)]

        Bfmathalfsq = [Bfmathalf[ff].dot(Bfmathalf[ff])
            for ff in xrange(punc)]

    # --------------
    # --------------
    # Get the Eigenvalues and vectors
    if onlyintegers:
        eigen_folder = 'results/eigendata/Np_'+ str(punc) + '_i_'+ Eigen_date +'/'
    else:
        eigen_folder = 'results/eigendata/Np_'+ str(punc) + '_' + Eigen_date +'/'

    filename = 'LT_' + str(lattice_type) + '_Trp_' + str(TOPRIGHT) + '_K_'+ str(Kspec) + 'Ktil_' + str(Ktilde) + 'vectors'

    Wavefuncs = np.load(eigen_folder+filename+'.npy')

    filename = 'LT_' + str(lattice_type) + '_Trp_' + str(TOPRIGHT) + '_K_'+ str(Kspec) + 'Ktil_' + str(Ktilde)

    Groundstate = np.load(eigen_folder+filename+'.npy')

    
    # Normalize to energy per puncture
    Groundstate /= float(punc)

    # How many energy levels
    TotalEigs = len(Groundstate[0])

    # How many intervals (number of points -1)
    maxin = len(Groundstate) - 1
    step = (alp_max-alp_min)/float(maxin)
    # alpha goes from alpmin to alpmx
    alps = np.arange(alp_min, alp_max + step, step)


    # --------------
    # --------------
    # Calculate the gap
    # Calculate the zeroth energy
    print '--------'
    print 'Gap calculation and derivatives'
    zeroenergy=[]
    zeroenergywvfn = []
    for nn, grounds in enumerate(Groundstate):
        posmin = np.where(grounds == grounds.min())[0][0]
        zeroenergy.append(grounds[posmin])
        # to calculate order parameter with ground state wave function
        zeroenergywvfn.append(Wavefuncs[nn][:,posmin])

    # Sort the Ground states
    Groundstatetemp=[]
    for grounds in Groundstate:
         Groundstatetemp.append(np.sort(grounds))

    Groundstate = Groundstatetemp

    # Calculate gap value
    energygap = [np.array([Groundstate[nn][mm] for nn in xrange(maxin + 1)])
                 - np.array(zeroenergy) for mm in xrange(TotalEigs)]


    # --------------
    # --------------
    # Create derivative operators
    firstder_oper = -0.5*np.eye(maxin+1,k=-1)+0.5*np.eye(maxin+1,k=1)
    firstder_oper /= step

    secondder_oper = -2*np.eye(maxin+1)+np.eye(maxin+1,k=1)+np.eye(maxin+1,k=-1)
    secondder_oper /= step**2

    # Calculate first derivative
    d1zero = np.dot(firstder_oper,zeroenergy)

    # Calculate second derivative
    d2zero = np.dot(secondder_oper,zeroenergy)

    critic = np.where(d2zero == d2zero[1:-1].min())[0]


    # --------------
    # --------------
    # Do the order parameter calculation
    print '--------'
    print 'Order parameter calculation'
    order =[]
    ordersq = []
    for nn, alp in enumerate(alps):
        magn = np.vdot(zeroenergywvfn[nn],zeroenergywvfn[nn])
        order.append(
            sum([
                np.vdot(zeroenergywvfn[nn],Bfmathalf[ff].dot(zeroenergywvfn[nn]))/magn
                for ff in xrange(punc)]))
        ordersq.append(
            sum([
                np.vdot(zeroenergywvfn[nn],Bfmathalfsq[ff].dot(zeroenergywvfn[nn]))/magn
                for ff in xrange(punc)
                ]))

    # The result of the expectation value is a 1dim matrix (which is supposedly real since W(1/2) is hermitian)
    order = np.array([np.real(order[nn]) for nn in xrange(maxin + 1)])/punc

    ordersq = np.array([np.real(ordersq[nn]) for nn in xrange(maxin + 1)])/punc
    
    
    # Fourth order Binder cumulant

    binder = 1 - (1./3.)* np.power(order,4) / np.power(ordersq,2)

    # Susceptibility

    Chi = punc*(ordersq - np.power(order,2)) / alps

    # Spec heat 

    # --------------
    # --------------
    # Locate the critical point
    print '--------'
    print 'Estimate point calculation'

    morethanone = ''
    for nn, ind in enumerate(critic):
        print 'Critical point at:', alps[ind]
        if nn > 0:
            print 'Warning: more of one minimum found'
            morethanone = '*'
        else:
            alpcrit = alps[ind]

    # Save the critical value and the energy gap

    # crit_file = 'results/critpoint.txt'

    # alreadywritten = ''

    # with open(crit_file, 'r') as f:
    #     for line in f:
    #         alreadywritten += line

    # with open(crit_file+'.tmp', 'w') as f:
        
    #     if onlyintegers: onin = 'T'
    #     else: onin = 'F'

    #     f.write(alreadywritten+'\n')

    #     f.write(str(punc) + ' '+ str(alpcrit)  + ' ' + str(lattice_type) + ' ' + onin + ' ' +
    #      str(Kspec)+ ' ' + str(Ktilde)+ ' '+ morethanone)

    # os.rename(crit_file+'.tmp',crit_file)


    print 'Saving the results to plot later'
    if not os.path.exists(save_folder):
        os.makedirs(save_folder)

    name0 = 'Torus'+str(punc)+str(lattice_type)+'_'+str(Kspec)+'_'+str(Ktilde)+'_alps'
    name1 = 'Torus'+str(punc)+str(lattice_type)+'_'+str(Kspec)+'_'+str(Ktilde)
    name2 = 'Torus_de1_'+str(punc)+str(lattice_type)+'_'+str(Kspec)+'_'+str(Ktilde)
    name3 = 'Torus_de2_'+str(punc)+str(lattice_type)+'_'+str(Kspec)+'_'+str(Ktilde)
    name4 = 'Torus_gap_'+str(punc)+str(lattice_type)+'_'+str(Kspec)+'_'+str(Ktilde)
    name5 = 'Torus_order_'+str(punc)+str(lattice_type)+'_'+str(Kspec)+'_'+str(Ktilde)
    name6 = 'Torus_binder_'+str(punc)+str(lattice_type)+'_'+str(Kspec)+'_'+str(Ktilde)
    name7 = 'Torus_susc_'+str(punc)+str(lattice_type)+'_'+str(Kspec)+'_'+str(Ktilde)

    np.save(save_folder+name0,alps)
    np.save(save_folder+name1,Groundstate)
    np.save(save_folder+name2,d1zero)
    np.save(save_folder+name3,[d2zero,alpcrit])
    np.save(save_folder+name4,energygap)
    np.save(save_folder+name5,order)
    np.save(save_folder+name6,binder)
    np.save(save_folder+name7,Chi)


if __name__ == '__main__':

    if len(sys.argv) > 2:
        raise ValueError('Invalid number of arguments')
    else:
        process()



