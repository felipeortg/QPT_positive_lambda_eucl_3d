#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2018-04-20 10:31:18
# @Author  : Felipe Ortega (felipeortegagama@gmail.com)
# @Title   : Create kin Hilbert

import os
from itertools import combinations as itcomb
from config import *
import quantum_defs as qdfs
from edgelabeler import edgelabeler_nonadj

# Translates from a number n from base 10 to any base b
# Stores it in a list of specific length
def numberToBase(n, b, lis_length):
    if n == 0:
        return [0 for nn in xrange(lis_length)]
    digits = []
    while n:
        digits.append(int(n % b))
        n /= b 
    return [0 for nn in xrange(lis_length-len(digits))] + digits[::-1]

#--------------
# Define the "admissible" function, allowed states in torus.
def Admiss(State, vertices, K):
    # Check the coupling rules
    if (qdfs.prod([qdfs.Delt(
            State[vertices[nn][0]],
            State[vertices[nn][1]],
            State[vertices[nn][2]],
            K)
            for nn in xrange(len(vertices)) ])):
        return True
    else:
        return False


# Optimized algorithm to calculate Hilbert space for K=1 (only closed string-nets)
def findcorrect(zz,vertices,edgecounter):
    tolvertices = len(vertices)
    perm = itcomb(range(edgecounter), zz)
    correctconf_fun = []
    for tupleconf in perm:
        correctvertex = 0
        settuple = set(tupleconf)
        for vv in vertices:
            if len(set(vv)&settuple)%2: 
                break
            correctvertex += 1
        if tolvertices == correctvertex:
            correctconf_fun.append(tupleconf)
    return correctconf_fun


def getTJk3(start, stop, vertices, edgecounter):
    tolvertices = len(vertices)
    correctconf_fun = []
    for nn in xrange(start, stop):

        correctvertex = 0
        check = numberToBase(nn, 2, edgecounter)
        for vv in vertices:
            if sum([check[aa] for aa in vv]) == 1:
                break
            correctvertex += 1
        if tolvertices == correctvertex:
            correctconf_fun.append(check)

    return 2 * np.array(correctconf_fun)


    # return [ 2 * np.array( numberToBase(nn, 2, edgecounter) ) 
    #         for nn in xrange(start, stop)
    #         if Admiss( 2 * np.array( numberToBase(nn, 2, edgecounter) ), vertices, Kspec)]


def getTJ(start, stop, vertices, edgecounter):
    return [np.array( numberToBase(nn, Kspec+1, edgecounter) ) 
            for nn in xrange(start, stop)
            if Admiss( np.array( numberToBase(nn, Kspec+1, edgecounter) ), vertices, Kspec)]


def generatehilbe():

    # Label top right is the number of the face in the top right of face 1
    # Generate the labeling of the lattice
    #--------

    vertices, InFace, OutFace, edgecounter = edgelabeler_nonadj(punc,TOPRIGHT)
    

    tolvertices = len(vertices)
    print 'vertices: ', tolvertices
    print 'Calculating the allowed Hilbert space'

    # Improved way to find Admiss states (still exponential time), only for K=1, K=2i
    if (Kspec == 1) or (Kspec == 2 and onlyintegers):

        cores = multiprocessing.cpu_count()
        pool = multiprocessing.Pool( cores )

        allowed = [zz for zz in xrange(0, 2*punc+1) if zz%2==0 ]

        results = [pool.apply_async( findcorrect, (zz,vertices,edgecounter) ) for zz in allowed]

        correctconfpar = [conf
                            for result in results
                            for conf in result.get()]

        TJ=[]
        for nn in correctconfpar:
            TJ.append([int(xx in nn) for xx in xrange(edgecounter)])

        if onlyintegers and Kspec == 2:
            TJ = [[2*spin for spin in conf] for conf in TJ]

    
    # PD: old code moved to the end

    elif Kspec==3 and onlyintegers:
        #--------------
        # K 3 only integers 
        print 'Only using states with integers'

        # Check only integer configurations:
        kinemsize = 2**edgecounter

        # TJ=[]
        # for nn in xrange(kinemsize):
        #     if nn%10000 == 0:
        #         print nn
        #     if Admiss( 2 * np.array( numberToBase(nn, 2, edgecounter) ), vertices, Kspec):
        #         TJ.append(2 * np.array( numberToBase(nn, 2, edgecounter) ))


        # Generate states as if kinemsize K=1 then just multiply the state by 2
        
        cores = multiprocessing.cpu_count()
        pool = multiprocessing.Pool( cores )

        TASK = [(int(nn*kinemsize/(cores*3.)), int((nn+1)*kinemsize/(cores*3.)),vertices, edgecounter) for nn in range(cores*3)]

        
        results = [pool.apply_async( getTJk3, t ) for t in TASK]

        TJ = [state
                for states in results
                for state in states.get()]

        # for nn in TJ:
        #     print nn

        # print len(TJ)
        # PD: old code moved to the end


    # The most general procedure to generate the Hilbert space
    else:

        kinemsize = (Kspec+1)**edgecounter

        cores = multiprocessing.cpu_count()
        pool = multiprocessing.Pool( cores )

        TASK = [(int(nn*kinemsize/(cores*3.)), int((nn+1)*kinemsize/(cores*3.)),vertices, edgecounter) for nn in range(cores*3)]
        
        results = [pool.apply_async( getTJ, t ) for t in TASK]

        TJ = [state
                for states in results
                for state in states.get()]



    

    # Save the Hilbert space

    print 'Save Hilbert vectors'

    if onlyintegers:
        save_folder = 'results/hilbe/Np_'+ str(punc) + '_i/'
    else:
        save_folder = 'results/hilbe/Np_'+ str(punc) + '/'

    if not os.path.exists(save_folder):
        os.makedirs(save_folder)


    filename = 'LT_' + str(lattice_type) + '_Trp_' + str(TOPRIGHT) + '_K_'+ str(Kspec)

    Hilbevectors = [TJ, vertices, InFace, OutFace, edgecounter]

    np.save(save_folder+filename,Hilbevectors) 

    return TJ, vertices, InFace, OutFace, edgecounter


if __name__ == '__main__':

    if len(sys.argv) > 1:
        raise ValueError('Invalid number of arguments')
    else:
        generatehilbe()



        # end = timer()

        # print(end- start)

        # start = timer() 

        # # Iterables with the position of the non zero edges
        # possiblecombin = []
        # for zz in xrange(0, 2*punc+1):
        #     if zz%2==0: # only an even number of colored edges is actually admissible
        #         possiblecombin.append(itcomb(range(edgecounter), zz))
 
        # correctconf = []

        # for mm in xrange(len(possiblecombin)):
        #     print 'perm', mm
        #     for tupleconf in possiblecombin[mm]:

        #         correctvertex = 0
                
        #         # Check that it satisfies Admiss at all vertices
        #         for vv in vertices:
        #             # for K=1 there can only be 2 or 0 spins per vertex
        #             # The admissible function simplifies in this condition
        #             if len(set(vv)&set(tupleconf))%2: 
        #                 break
        #             correctvertex += 1

        #         if tolvertices == correctvertex:
        #             correctconf.append(tupleconf)

        # end = timer()

        # print(end- start)        

        # # Fill up the TJ list with the correct configurations
        # TJ=[]
        # for nn in correctconf:
        #     TJ.append([int(xx in nn) for xx in xrange(edgecounter)])

        # print 'nonpar', len(TJ)
        # print 'par', len(TJpar)




        # TJ = [ 2 * np.array( numberToBase(nn, 2, edgecounter) ) 
        #     for nn in xrange(kinemsize)
        #     if Admiss( 2 * np.array( numberToBase(nn, 2, edgecounter) ), vertices, Kspec)]

        # end = timer()

        # print(end- start)

        # print len(TJ)
        # print len(TJpar)




        # TJ = [np.array( numberToBase(nn, Kspec+1, edgecounter) ) 
        #     for nn in xrange(kinemsize)
        #     if Admiss( np.array( numberToBase(nn, Kspec+1, edgecounter) ), vertices, Kspec)]




