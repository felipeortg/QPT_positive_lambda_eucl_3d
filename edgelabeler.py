#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2018-04-11 16:01:38
# @Author  : Felipe Ortega (felipeortegagama@gmail.com)
# @Title   : Edgelabeler

# Label top right is the number of the face in the top right of face 1

# Generate the labeling of the lattice
#--------

from config import *


def edgelabeler_nonadj(Numberpunctures,TOPRIGHT):

    if lattice_type == 'sq':
        if punc == 4:
            # Generate the labeling of the lattice
            #--------
            InFace = [
                [0,1,2,3,4,5],
                [3,6,7,0,8,9],
                [10,4,9,11,1,7],
                [11,8,5,10,6,2]]

            OutFace = [
                [7,11,6,9,10,8],
                [2,10,1,5,11,4],
                [5,3,8,2,0,6],
                [9,0,4,7,3,1]]

            vertices = [
            [0,1,7],
            [1,2,11],
            [2,3,6],
            [3,4,9],
            [4,5,10],
            [5,0,8],
            [6,7,10],
            [8,9,11]]

            edgecounter = 12

        elif punc == 9:

            # Generate the labeling of the lattice
            #--------
            InFace = [
                [0,1,2,3,4,5],
                [3,6,7,8,9,10],
                [8,11,12,0,13,14],
                [15,16,17,18,6,2],
                [18,19,20,21,11,7],
                [21,22,23,15,1,12],
                [24,13,5,25,19,17],
                [25,4,10,26,22,20],
                [26,9,14,24,16,23]]

            OutFace = [
                [12,15,6,10,25,13],
                [2,18,11,14,26,4],
                [7,21,1,5,24,9],
                [23,24,19,7,3,1],
                [17,25,22,12,8,6],
                [20,26,16,2,0,11],
                [14,0,4,20,18,16],
                [5,3,9,23,21,19],
                [10,8,13,17,15,22]]

            vertices = [
            [0,1,12],
            [1,2,15],
            [2,3,6],
            [3,4,10],
            [4,5,25],
            [5,13,0],
            [6,7,18],
            [7,8,11],
            [8,9,14],
            [9,10,26],
            [11,12,21],
            [13,14,24],
            [15,16,23],
            [16,17,24],
            [17,18,19],
            [19,20,25],
            [20,21,22],
            [22,23,26]]

            edgecounter = 27

        else:
            print 'Error: Only square type lattices are with 4 and 9 punctures'

            sys.exit()

    else:

        numfaces,topright = Numberpunctures, TOPRIGHT
        topleft = numfaces-topright + 1

        if (topright==1 or topright==numfaces or 
            topleft==1 or topleft==numfaces):# redundant, last one implies this
            raise ValueError('No lattices with autoadjacent faces')

        ## See diagram to understand ma1, ma, nm1a, nma
        # Association of Js and Ls per face (i.e. the face 1 is surrounded by the edges 0, 1 etc; and the edges "outside" of the face 1 are ...).
        InFace = []
        edgecounter = 0 # first edge is '0'
        for ff in xrange(numfaces):
            InFace.append([])
            thisface = InFace[ff]

            #top edge
            if ff == 0:
                thisface.append(edgecounter)
                edgecounter += 1
            else:
                thisface.append(InFace[ff-1][3])

            # top right edge
            ma1 = (topright+ff-1)%numfaces
            if ma1 > ff-1: # condition is ma1 bigger or equal than ff
                thisface.append(edgecounter)
                edgecounter += 1

            else:
                thisface.append(InFace[ma1][4])

            # bot right edge
            ma = (ma1+1)%numfaces
            if ma > ff-1: # condition is ma1 bigger or equal than ff
                thisface.append(edgecounter)
                edgecounter += 1
            else:
                thisface.append(InFace[ma][5])      

            # bot edge
            if (ff+1)%numfaces > ff:
                thisface.append(edgecounter)
                edgecounter += 1
            else:
                thisface.append(InFace[0][0])  

            # bot left edge
            nm1a = (topleft + ff)%numfaces # labels begin with 0, human numbering differs by 1
            if nm1a > ff: # the condition is only bigger since the equal has been labeled
                thisface.append(edgecounter)
                edgecounter += 1
            else:
                thisface.append(InFace[nm1a][1])

            # top left edge
            nma = (nm1a - 1)%numfaces
            if nma > ff: # the condition is only bigger since the equal has been labeled
                thisface.append(edgecounter)
                edgecounter += 1
            else:
                thisface.append(InFace[nma][2])

        # See diagram to understand ma1, ma, nm1a, nma
        OutFace = []
        for ff in xrange(numfaces):
            
            ma1 = (topright+ff-1)%numfaces # top right face
            ma = (ma1+1)%numfaces # bot right face
            nm1a = (topleft+ff)%numfaces # bot left face
            nma = (nm1a - 1)%numfaces # top left face

            OutFace.append(
                [InFace[ma1][5],
                InFace[ma1][3],
                InFace[ma][4],
                InFace[nm1a][2],
                InFace[nma][3],
                InFace[nma][1]
                ])

        vertices = []
        countedvert = [[True for nn in xrange(6)] for mm in xrange(numfaces)]

        for ff in xrange(numfaces):

            ma1 = (topright+ff-1)%numfaces # top right face
            ma = (ma1+1)%numfaces # bot right face
            nm1a = (topleft+ff)%numfaces # bot left face
            nma = (nm1a - 1)%numfaces # top left face
            adjfaces = [(ff-1)%numfaces,ma1,ma,(ff+1)%numfaces,nm1a,nma]
            for vert in xrange(6):
                if countedvert[ff][vert]:
                    vertices.append([
                        InFace[ff][vert],
                        InFace[ff][(vert+1)%6],
                        OutFace[ff][vert]])
                    countedvert[ff][vert] = False
                    countedvert[adjfaces[vert]][(vert+2)%6] = False
                    countedvert[adjfaces[(vert+1)%6]][(vert+4)%6] = False

    return vertices, InFace, OutFace, edgecounter

#Vectors generated: "vertices, InFace, OutFace", the number of edges: "edgecounter"

