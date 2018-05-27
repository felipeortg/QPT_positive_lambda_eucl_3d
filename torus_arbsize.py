# -*- coding: utf-8 -*-
# @Date    : 2018-03-27 16:10:24
# @Author  : Felipe Ortega (felipeortegagama@gmail.com)
# @Title   : torus


import operator
import quantum_defs as qdfs
import gen_hil as ghil
from config import *

# sys.exit()


#--------------
# SPARSE CALCULATION DEFINITION
# Needed to be here to parallelize


# Find the invariant subspace with same config out of the face
# Receives the Hilbert space vectors configs as input
# The list with the outedges for the given face is also needed
def findsubspaces(face, HVs, outedges_ff):
    NT = len(HVs) # Hilbert space size
    sub = [[ 0 ]] # array to save subspaces of ff
    leadst = [0] # states to compare with (1 per subspace)
    leadstout = [[HVs[0][nn] for nn in outedges_ff]]
    
    for aa in xrange(1, NT): # classify all states (the first state is done)
        classified = False
        outtoclass = [HVs[aa][nn] for nn in outedges_ff]
        
        for mm, ls in enumerate(leadst): # check if it matches a subspace
            
            if outtoclass == leadstout[mm]:

                sub[mm].append(aa)
                
                classified = True
                
                break

        if not classified:
            sub.append([ aa ])
            leadst.append(aa)
            leadstout.append(outtoclass)

    return sub


# Flatness operator definition
# B-matrix (per subspace):
def Bfacesparalel(subs, Face, outedges, InFace, OutFace, K):
    # method checks for out-edges before calculating ME
    Edges = 6  # ATTENTION: just works with 6-faced hexagon
    return  [ 
                (qdfs.DT(K)**(-2) * 
                sum([
                    qdfs.vi(mm, K)**2 
                    * qdfs.prod([
                        qdfs.FS(
                            StateN[ InFace[Face][(nn+1) % Edges] ], 
                            StateN[ OutFace[Face][nn] ], 
                            StateN[ InFace[Face][nn] ],
                            StateJ[ InFace[Face][nn] ], 
                            mm, 
                            StateJ[ InFace[Face][(nn+1) % Edges] ],
                            K
                            )
                        for nn in xrange(Edges)
                        ])
                    for mm in xrange(K+1)
                    ]))
                    for StateN in subs
                    for StateJ in subs]


# Order parameter operator definition
# B-matrix-half (per subspace):
def Bface_half_paralel(subs, Face, outedges, InFace, OutFace, K):
    if onlyintegers:
        m = 2 # Notation of twice the spin to work only with integers
    else:
        m = 1 # Notation of twice the spin to work only with integers

    # method checks for out-edges before calculating ME
    Edges = 6  # ATTENTION: just works with 6-faced hexagon
    return  [(1/(qdfs.vi(m, K)**2)
                * qdfs.prod([
                    qdfs.FS(
                        StateN[ InFace[Face][(nn+1) % Edges] ], 
                        StateN[ OutFace[Face][nn] ], 
                        StateN[ InFace[Face][nn] ],
                        StateJ[ InFace[Face][nn] ], 
                        m, 
                        StateJ[ InFace[Face][(nn+1) % Edges] ],
                        K
                        )
                    for nn in xrange(Edges)
                    ]))
                for StateN in subs
                for StateJ in subs]


# Input the state numbers    
def Bfacetildes(statenumN, statenumJ, Face, K):
    # method checks for out-edges before calculating ME
    if Admissmemo(statenumN) and Admissmemo(statenumJ):
        StateN, StateJ = TJ[statenumN], TJ[statenumJ]  
        return (qdfs.DT(K)**(-2) * 
        sum([
            qdfs.vi(mm, K)**2 
            * qdfs.prod([
                qdfs.FS(
                    StateN[ InFace[Face][(nn+1) % Edges] ],
                    StateN[ OutFace[Face][nn] ], 
                    StateN[ InFace[Face][nn] ],
                    StateJ[ InFace[Face][nn] ], 
                    mm, 
                    StateJ[ InFace[Face][(nn+1) % Edges ] ],
                    K
                    )
                for nn in xrange(Edges)
                ])
            for mm in xrange(K+1)
            ]))

    else:
        return 0



def torus_arbpuncspectrum():

    # Print parameters:

    print 'Torus with', punc, 'punctures/faces'
    print 'Lattice type', lattice_type
    print 'K =', Kspec,' and Ktilde =', Ktilde

    # Functions for the arb size torus (when no faces are autoadjacent)
    # -----------------

    # Label top right is the number of the face in the top right of face 1
    # Generate the labeling of the lattice
    # Get the Hilbert space vectors

    if onlyintegers:
        hilbe_folder = 'results/hilbe/Np_'+ str(punc) + '_i/'
    else:
        hilbe_folder = 'results/hilbe/Np_'+ str(punc) + '/'

    filename = 'LT_' + str(lattice_type) + '_Trp_' + str(TOPRIGHT) + '_K_'+ str(Kspec) +'.npy'

    # sys.exit()

    if not os.path.isfile(hilbe_folder+ filename): # check if the Hilbert space has been already calculated
        TJ, vertices, InFace, OutFace, edgecounter = ghil.generatehilbe()

    else: # get the Hilbe vectors
        print 'Hilbert vectors in memory, skipping the search of H_kin'
        Hilbevectors = np.load(hilbe_folder+filename)

        TJ = Hilbevectors[0]
        vertices  = Hilbevectors[1]
        InFace = Hilbevectors[2]
        OutFace = Hilbevectors[3]
        edgecounter = Hilbevectors[4]

    #Vectors generated: "vertices, InFace, OutFace", the number of edges: "edgecounter"
   
    # vector total edges to verify later if ME are 0 just by the spins out-of-face
    TotalEdges = range(edgecounter)
    print 'Number of edges:', edgecounter

    # --------------
    # --------------
    # Admissible for Bp tilde
    # Memoized: define the "admissible" function, allowed states in torus, useful for k-tilde
    
    if Ktilde == 0:
        @qdfs.memoize_sing # only going to use it with Ktilde = 0 anyway, and memoize_sing is faster
        def Admissmemo(statenum):
            if statenum == 0:
                return True
            else:
                return False

    else:
        @qdfs.memoize_sing # only going to use it with Ktilde anyway, and memoize_sing is faster
        def Admissmemo(statenum):

            State = TJ[statenum]
            if (qdfs.prod([qdfs.Delt(
                    State[vertices[nn][0]],
                    State[vertices[nn][1]],
                    State[vertices[nn][2]],
                    Ktilde)
                    for nn in xrange(len(vertices)) ])):
                return True
            else:
                return False

    # --------------
    # --------------
    # Print Hilbert space to check it

    NT = len(TJ)

    for nn in xrange(NT):
        print np.array(TJ[nn])/2.
        
        if nn > 15:
            break # don't overload the log file

    print 'Dimension of H :', NT

    

    # --------------
    # --------------
    # Check if the ME are already calculated
    if onlyintegers:
        matrix_folder = 'results/matrices/Np_'+ str(punc) + '_i_'+ B_date + '/'
    else:
        matrix_folder = 'results/matrices/Np_'+ str(punc) + '_'+ B_date  + '/'
    
    filename = 'LT_' + str(lattice_type) + '_Trp_' + str(TOPRIGHT) + '_K_'+ str(Kspec) + '.npy'

    # Flag to know if Bp(Kspec) is in memory
    bp_k_mat = os.path.isfile(matrix_folder+ filename)

    filename = filename[:-4] + '_Ktil_' + str(Ktilde) +'_tBp.npy'
    
    # Flag to know if Bp(Ktilde) is in memory
    tilbp_tilk_mat = os.path.isfile(matrix_folder+ filename)

    if bp_k_mat and tilbp_tilk_mat:

        print 'Matrix elements in memory, skipping their explicit calculation'

        sys.exit()

    elif bp_k_mat:

        print 'Matrix tilde missing, doing its calculation'


    # Faces in a torus are hexagons and have 6 edges
    Edges = 6
    Faces = punc

    # which spins are outside of each plaquette 
    # (Bp acts on the subspace of same spins out of plaquette)
    outedges = [list(set(TotalEdges)-set(InFace[nn])) for nn in xrange(Faces)]

    if NT < 2**10+1: # 'dense exact diagonalization procedure'
        print 'Dense calculation'

        #--------------
        # Flatness operator definition

        # B-matrix element:
        def Bface(StateN, StateJ, Face, K):
            if [StateN[nn] for nn in outedges[Face]] == [StateJ[nn] for nn in outedges[Face]]:
                return (qdfs.DT(K)**(-2) * 
                sum([
                    qdfs.vi(mm, K)**2 
                    * qdfs.prod([
                        qdfs.FS(
                            StateN[ InFace[Face][(nn+1) % Edges] ],
                            StateN[ OutFace[Face][nn] ], 
                            StateN[ InFace[Face][nn] ],
                            StateJ[ InFace[Face][nn] ], 
                            mm, 
                            StateJ[ InFace[Face][(nn+1) % Edges ] ],
                            K
                            )
                        for nn in xrange(Edges)
                        ])
                    for mm in xrange(K+1)
                    ]))
            else:
                return 0

        # Input the state numbers
        def Bfacetilde(statenumN, statenumJ, Face, K):
            if Admissmemo(statenumN) and Admissmemo(statenumJ):
                StateN, StateJ = TJ[statenumN], TJ[statenumJ]
                if [StateN[nn] for nn in outedges[Face]] == [StateJ[nn] for nn in outedges[Face]]:
                    return (qdfs.DT(K)**(-2) * 
                    sum([
                        qdfs.vi(mm, K)**2 
                        * qdfs.prod([
                            qdfs.FS(
                                StateN[ InFace[Face][(nn+1) % Edges] ],
                                StateN[ OutFace[Face][nn] ], 
                                StateN[ InFace[Face][nn] ],
                                StateJ[ InFace[Face][nn] ], 
                                mm, 
                                StateJ[ InFace[Face][(nn+1) % Edges ] ],
                                K
                                )
                            for nn in xrange(Edges)
                            ])
                        for mm in xrange(K+1)
                        ]))
                else:
                    return 0
            else:
                return 0

        # Order parameter operator definition
        # B-matrix-half (per subspace):
        def Bface_half(StateN, StateJ, Face, K):
            if onlyintegers:
                m = 2 # Notation of twice the spin to work only with integers
            else:
                m = 1 # Notation of twice the spin to work only with integers

            if [StateN[nn] for nn in outedges[Face]] == [StateJ[nn] for nn in outedges[Face]]:
                return (1/(qdfs.vi(m, K)**2)
                    * qdfs.prod([
                        qdfs.FS(
                            StateN[ InFace[Face][(nn+1) % Edges] ],
                            StateN[ OutFace[Face][nn] ], 
                            StateN[ InFace[Face][nn] ],
                            StateJ[ InFace[Face][nn] ], 
                            m,
                            StateJ[ InFace[Face][(nn+1) % Edges ] ],
                            K
                            )
                        for nn in xrange(Edges)
                        ]))
            else:
                return 0


        #--------------
        #--------------
        # Calculate Hamiltonian and print
        print 'Begin ME explicit calculation'



        # Bfmat = []
        
        # for ff in xrange(Faces):
        #     temp1 =[]
        #     for mm in xrange(NT):
        #         temp2 =[]
        #         for nn in xrange(NT):

        #             temp2.append(Bface(TJ[nn],TJ[mm],ff,Kspec))

        #         temp1.append(temp2)

        #     Bfmat.append(temp1)

        if not bp_k_mat:
            Bfmat = [
                np.array([[Bface(TJ[nn],TJ[mm],ff,Kspec) for nn in xrange(NT)] 
                    for mm in xrange(NT)])
            for ff in xrange(Faces)] 

            print 'Bp done'

            Bfmathalf = [
                np.array([[Bface_half(TJ[nn],TJ[mm],ff,Kspec) for nn in xrange(NT)] 
                    for mm in xrange(NT)])
            for ff in xrange(Faces)]

            print 'B-half done' 

        Bfmattilde = [
            np.array([[Bfacetilde(nn,mm,ff,Ktilde) for nn in xrange(NT)] 
                for mm in xrange(NT)])
        for ff in xrange(Faces)] 

        print 'Bptilde done, saving the values'        

        #--------------
        #--------------

        # Save matrices
        if not bp_k_mat:

            if onlyintegers:
                save_folder = 'results/matrices/Np_'+ str(punc) + '_i_'+str(dt.datetime.now())[:-16] +'/'
            else:
                save_folder = 'results/matrices/Np_'+ str(punc) + '_' +str(dt.datetime.now())[:-16] +'/'

            if not os.path.exists(save_folder):
                os.makedirs(save_folder)

            filename = 'LT_' + str(lattice_type) + '_Trp_' + str(TOPRIGHT) + '_K_'+ str(Kspec) 

            np.save(save_folder+filename, ['dense'])

            filename = filename + '_Bp'

            np.save(save_folder+filename,Bfmat)

            filename = filename[:-2] + 'B-half'

            np.save(save_folder+filename,Bfmathalf)

            filename = filename[:-6] + 'Ktil_' + str(Ktilde) +'_tBp'
        else:

            if onlyintegers:
                save_folder = 'results/matrices/Np_'+ str(punc) + '_i_'+B_date +'/'
            else:
                save_folder = 'results/matrices/Np_'+ str(punc) + '_' +B_date +'/'
            
            filename = 'LT_' + str(lattice_type) + '_Trp_' + str(TOPRIGHT) + '_K_'+ str(Kspec)
            filename += '_Ktil_' + str(Ktilde) +'_tBp'

        np.save(save_folder+filename,Bfmattilde)

    # ----------
    # ----------
    # ----------
    # Sparse matrix methods
    else:
        print 'Sparse calculation'

        #--------------
        # Flatness operator definition (the other one is at the beginning)
        # Input the state numbers    
        def Bfacetildes(statenumN, statenumJ, Face, K):
            # method checks for out-edges before calculating ME
            if Admissmemo(statenumN) and Admissmemo(statenumJ):
                StateN, StateJ = TJ[statenumN], TJ[statenumJ]  
                return (qdfs.DT(K)**(-2) * 
                sum([
                    qdfs.vi(mm, K)**2 
                    * qdfs.prod([
                        qdfs.FS(
                            StateN[ InFace[Face][(nn+1) % Edges] ],
                            StateN[ OutFace[Face][nn] ], 
                            StateN[ InFace[Face][nn] ],
                            StateJ[ InFace[Face][nn] ], 
                            mm, 
                            StateJ[ InFace[Face][(nn+1) % Edges ] ],
                            K
                            )
                        for nn in xrange(Edges)
                        ])
                    for mm in xrange(K+1)
                    ]))

            else:
                return 0

        # Find the subspaces of the Bp operators
        print 'Find subspaces of Bp operators'

        # --------------
        # Multiprocessing definitions

        cores = multiprocessing.cpu_count()
        pool = multiprocessing.Pool( cores )

        # Create the workers
        results = [
            pool.apply_async( findsubspaces , (ff, TJ, outedges[ff]))
        for ff in xrange(Faces)]

        # Get the subspaces

        subspaces = [subs.get() for subs in results]
        

        # Input for parallelization is easier as configurations (I think since the TJ wont be shared)
        subspacesTJ = [[[TJ[nn]
                for nn in subsp]
            for subsp in subspaces[ff]]
        for ff in xrange(Faces)]

        #--------------
        #--------------
        # Calculate Hamiltonian and print
        print 'Begin ME explicit calculation'

        if not bp_k_mat:

            results = [[
                        pool.apply_async( Bfacesparalel, (subs, ff, outedges, InFace, OutFace, Kspec))
                        for subs in subspacesTJ[ff]]
                    for ff in xrange(Faces)]

            
            Bfmatdata = [[ME_of_Bf

                        for suboff in faceresult
                        for ME_of_Bf in suboff.get()]

                        for faceresult in results]


            print 'Bp done'


            results = [[
                        pool.apply_async( Bface_half_paralel, (subs, ff, outedges, InFace, OutFace, Kspec))
                        for subs in subspacesTJ[ff]]
                    for ff in xrange(Faces)]

            
            Bfmathalfdata = [[ME_of_Bfhalf

                        for suboff in faceresult
                        for ME_of_Bfhalf in suboff.get()]

                        for faceresult in results]

            print 'B-half done'



        Bfmattildata = [[
                    Bfacetildes(nn,mm,ff,Ktilde) # uses statenum to ease the Admissmemo
                    
                    for subs in subspaces[ff] 
                    for nn in subs
                    for mm in subs]

            for ff in xrange(Faces)]

        rowssparse = [[
                    nn 

                    for subs in subspaces[ff] 
                    for nn in subs
                    for mm in subs]

            for ff in xrange(Faces)]

        colssparse = [[
                    mm

                    for subs in subspaces[ff] 
                    for nn in subs
                    for mm in subs]

            for ff in xrange(Faces)]


        print 'Bptilde done, saving the values' 
        
        #--------------
        #--------------

        # Save matrices

        if not bp_k_mat:

            if onlyintegers:
                save_folder = 'results/matrices/Np_'+ str(punc) + '_i_'+str(dt.datetime.now())[:-16] +'/'
            else:
                save_folder = 'results/matrices/Np_'+ str(punc) + '_' +str(dt.datetime.now())[:-16] +'/'

            if not os.path.exists(save_folder):
                os.makedirs(save_folder)

            filename = 'LT_' + str(lattice_type) + '_Trp_' + str(TOPRIGHT) + '_K_'+ str(Kspec) 

            np.save(save_folder+filename, ['sparse'])

            filename = filename + '_Bp'

            np.save(save_folder+filename,[Bfmatdata,rowssparse,colssparse,NT])

            filename = filename[:-2] + 'B-half'

            np.save(save_folder+filename,[Bfmathalfdata,rowssparse,colssparse,NT])

            filename = filename[:-6] + 'Ktil_' + str(Ktilde) +'_tBp'

        # Only K-tilde calculation (save it in the same folder)
        else:

            if onlyintegers:
                save_folder = 'results/matrices/Np_'+ str(punc) + '_i_'+B_date +'/'
            else:
                save_folder = 'results/matrices/Np_'+ str(punc) + '_' +B_date +'/'

            # Avoid a too long line of code:    
            filename = 'LT_' + str(lattice_type) + '_Trp_' + str(TOPRIGHT) + '_K_'+ str(Kspec)
            filename += '_Ktil_' + str(Ktilde) +'_tBp'

        np.save(save_folder+filename,Bfmattildata)



if __name__ == '__main__':

    if len(sys.argv) > 1:
        raise ValueError('Invalid number of arguments')
    else:
        torus_arbpuncspectrum()


