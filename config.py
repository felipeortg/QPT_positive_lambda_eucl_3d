#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2018-04-19 19:20:15
# @Author  : Felipe Ortega (felipeortegagama@gmail.com)
# @Title   : Parameters

import os
import sys
import numpy as np
import datetime as dt
import multiprocessing
from timeit import default_timer as timer

# The lattice type will have the following convention
 
# For "linear" lattices letters 'a', 'b', 'c', ...
#   this letter depends on the label of the top right puncture of 1
#   'a': for top right puncture 1, i.e. lattice with adjacencies
#   'b', 'c', ...: for top right 2, 3, ...
 
# For squared matrices use 'sq'
 
# Lattice specifics
punc = int(10.0)
lattice_type = 'b'
TOPRIGHT = int(2.0)

# Hamiltonian and Hilbert space
onlyintegers = bool(False)
Kspec = int(1.0)
Ktilde = int(0.0)
if Ktilde != 0:
    print '------------'
    print '------------'
    print 'Warning: lattice dependent spectrum'
    print '------------'
    print '------------'

# Range of eigenvalue calculation
range_alp_def = int(1.0)

if range_alp_def == 0:
    alp_min = 0.0
    alp_max = 1.0
else:
    alp_min = 0.65
    alp_max = 0.74

# Number of calculated points
points_to_calc = int(225.0)

# Dates if only processing/plotting is needed
process_date = str(dt.datetime.now())[:-16] 
Eigen_date = '2018-05-16' # str(dt.datetime.now())[:-16] # '2018-04-19'
B_date = '2018-05-16' # str(dt.datetime.now())[:-16]     # '2018-04-19'


# Save (0) or plot (1)
save_show = int(0.0)
