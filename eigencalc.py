#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2018-04-19 15:50:59
# @Author  : Felipe Ortega (felipeortegagama@gmail.com)
# @Title   : Eigenvalues calculation


import scipy.sparse as spar
import scipy.sparse.linalg as spalg
from config import *


def parallel_eigensys(alpsublis, H1, H2, MaxEigs):
    # Avoid terminate the calculation at error
    APerror = False
    num_of_errors = 0
    consec_errors = 0
    errorat = []
    Groundstate = []
    Wavefuncs = []

    lensub = len(alpsublis)
    for ii, alp in enumerate(alpsublis):
        Hamiltonian = - alp*H1 - (1-alp)*H2
        #Careful with AP error
        try:
            eigenvalues, eigenvectors = spalg.eigsh(Hamiltonian, k = MaxEigs, which = 'LM', return_eigenvectors=True)

            # Debugging: raise an error at ii=2
            # if ii==2:
            #     raise spalg.eigen.ArpackError(1)

            if APerror == True and not (ii == 1): # won't work if APerror was found in ii = 0
                if consec_errors > 1:
                    raise RuntimeError('Too many Arpack errors')

                # ATTENTION:
                # Correction of Arpack errors done after full eigcalc 
                # (following commented lines are now useless)
                # Use last evaluation and current to replace missing values with mean
                # lastevals = (Groundstate[-1] + evals_large)/2.0
                # Groundstate.append(lastevals)

                # In the meantime extend twice current evals
                Groundstate.extend([eigenvalues, eigenvalues])
                Wavefuncs.extend([eigenvectors, eigenvectors])

                # change the state of the flags
                APerror = False
                consec_errors = 0

            elif APerror == True and (ii == 1):
                # if error happened on the first value extend twice the same value
                Groundstate.extend([eigenvalues, eigenvalues])
                Wavefuncs.extend([eigenvectors, eigenvectors])

                # change the state of the flags
                APerror = False
                consec_errors = 0
            
            else:
                Groundstate.append(eigenvalues)
                Wavefuncs.append(eigenvectors)
            
        except spalg.eigen.ArpackError as e:
            print e
            APerror = True
            errorat.append(alp)
            num_of_errors += 1
            consec_errors += 1

    # if last calculation was wrong re-append the last value
    if APerror == True:
        if consec_errors > 1:
            raise RuntimeError('Too many Arpack errors')

        Groundstate.append(Groundstate[-1])
        Wavefuncs.append(Wavefuncs[-1])

    
    # For the sake of impatience (and to make sure is doing stuff)
    print 'Subset done'

    return Groundstate, Wavefuncs, errorat


def eigenvaluescalc(points):

    # Get the matrices

    print 'Getting the ME from files'

    if onlyintegers:
        matrixfolder = 'results/matrices/Np_'+ str(punc) + '_i_'+ B_date +'/'
    else:
        matrixfolder = 'results/matrices/Np_'+ str(punc) + '_' + B_date +'/'

    # Get the Bps
    filename = 'LT_' + str(lattice_type) + '_Trp_' + str(TOPRIGHT) + '_K_'+ str(Kspec)

    mattype = np.load(matrixfolder+filename+'.npy')[0]

    if mattype == 'dense':

        filename +=  '_Bp'

        Bfmat = np.load(matrixfolder+filename+'.npy')

        # Size of the Hilbert space
        NT = len(Bfmat[0])

        # Get the Bps tilde
        filename = filename[:-2] + 'Ktil_' + str(Ktilde) +'_tBp'

        Bfmattilde = np.load(matrixfolder+filename+'.npy')

    elif mattype == 'sparse':

        filename +=  '_Bp'

        # Size of the Hilbert space
        Bfmatdata,rowssparse,colssparse,NT = np.load(matrixfolder+filename+'.npy')

        # Get the Bps tilde
        filename = filename[:-2] + 'Ktil_' + str(Ktilde) +'_tBp'

        Bfmattildata = np.load(matrixfolder+filename+'.npy')
        
        Bfmat = [
            spar.coo_matrix(
                (Bfmatdata[ff],
                    (rowssparse[ff],
                    colssparse[ff])
                ),
            shape = (NT,NT)).tocsr()

            for ff in xrange(punc)]

        Bfmattilde = [
            spar.coo_matrix(
                (Bfmattildata[ff],
                    (rowssparse[ff],
                    colssparse[ff])
                ),
            shape = (NT,NT)).tocsr()

            for ff in xrange(punc)]       


    # Do the eigenvalue calculation

    H1 = sum(Bfmat)
    H2 = sum(Bfmattilde)

    maxin = int(points)
    step = (alp_max- alp_min)/float(points)
    MaxEigs = 32 # for the case of sparse matrices

    # alpha goes from 0 to 1
    alps = np.arange(alp_min, alp_max + step, step)
    
    print 'Spectrum calculation (', points, '+ 1 points)'
    
    # Prepare for parallel computation

    cores = multiprocessing.cpu_count()
    pool = multiprocessing.Pool( cores )

    # Calculate eigenvalues
    if NT < 2**10+1:
        Groundstate = []
        Wavefuncs = []
        for nn, alp in enumerate(alps):
            Hamiltonian = - alp*H1 - (1-alp)*H2
            eigenvalues, eigenvectors = np.linalg.eigh(Hamiltonian)
            Groundstate.append(eigenvalues)
            Wavefuncs.append(eigenvectors)
            # For the sake of impatience (and to make sure is doing stuff)
            if nn%10 == 0:
                print nn

        num_of_errors = 0

    else:
        # Divide the work along the cores
        alpsize = maxin + 1

        # Divide the work in sets of 20 points
        alpsublis = [(alps[20*nn:20*nn+20], H1, H2, MaxEigs) for nn in range(alpsize/20)]

        # Add the last batch in case there are alps left
        if alpsize%20 != 0:
            alpsublis.append((alps[((alpsize/20) * 20):], H1, H2, MaxEigs))

        # Do the diagonalization
        results = [pool.apply_async( parallel_eigensys, aa ) for aa in alpsublis]

        Grst_err = [rr.get() for rr in results]

        errors = [err
                    for list_err in Grst_err
                    for err in list_err[2]
                    ]
        num_of_errors = len(errors)

        Wavefuncs = [wfcs 
                        for list_wfcs in Grst_err
                        for wfcs in list_wfcs[1]
                        ]

        Groundstate = [grst 
                        for list_grst in Grst_err
                        for grst in list_grst[0]
                        ]

        # Correction for ARPACK errors:
        print 'Corrections at :', errors
        for error in errors:# Corrections in case of Arpack errors
            # Will raise an IndexError in case of not finding the error in alps
            ind = np.where(alps == error)[0][0] # for the tuple and then the array 

            try:
                if ind == 0:
                    raise IndexError
                else:
                    Groundstate[ind] = (Groundstate[ind-1] + Groundstate[ind+1])/2.0
                    # The eigenvectors are not normalized anyway:
                    Wavefuncs[ind] = (Wavefuncs[ind-1] + Wavefuncs[ind+1])/2.0 

            except IndexError:
                print 'Warning: Error in first or last value (not corrected)'   

    # Save the Eigenvalues

    if onlyintegers:
        save_folder = 'results/eigendata/Np_'+ str(punc) + '_i_'+str(dt.datetime.now())[:-16] +'/'
    else:
        save_folder = 'results/eigendata/Np_'+ str(punc) + '_' +str(dt.datetime.now())[:-16] +'/'

    if not os.path.exists(save_folder):
        os.makedirs(save_folder)

    filename = 'LT_' + str(lattice_type) + '_Trp_' + str(TOPRIGHT) + '_K_'+ str(Kspec) + 'Ktil_' + str(Ktilde) + 'vectors'

    np.save(save_folder+filename,Wavefuncs)

    filename = 'LT_' + str(lattice_type) + '_Trp_' + str(TOPRIGHT) + '_K_'+ str(Kspec) + 'Ktil_' + str(Ktilde)

    np.save(save_folder+filename,Groundstate)

    if num_of_errors > 0:

        filename += '_errors' 
        np.save(save_folder+filename,errors)




if __name__ == '__main__':
    if len(sys.argv) > 3:
        raise ValueError('Invalid number of arguments')
    else:
        eigenvaluescalc(points_to_calc)



