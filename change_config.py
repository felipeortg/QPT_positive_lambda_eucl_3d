#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2018-04-19 21:45:55
# @Author  : Felipe Ortega (felipeortegagama@gmail.com)
# @Title   : Change config_file

import os
import sys

# -----------------
# Config file:
config_file = 'config.py'
counter = float(sys.argv[1])
lt_type = sys.argv[2]
top_right = float(sys.argv[3])
integers = bool(float(sys.argv[4])) # bool by itself recognizes '0' as True
kvalue = float(sys.argv[5])
ktild = float(sys.argv[6])
defalp = float(sys.argv[7])
points = float(sys.argv[8])
svnotplot = float(sys.argv[9])



text = ''

with open(config_file, 'r') as f:
    for line in f:
        if line[0:5] == 'punc ':
            text += 'punc = int(' + str(counter) +')\n'
            print 'Punctures: ' + str(counter)

        elif line[0:13] == 'lattice_type ':
            text += "lattice_type = '" + lt_type +"'\n"

        elif line[0:9] == 'TOPRIGHT ':
            text += 'TOPRIGHT = int(' + str(top_right) +')\n'

        elif line[0:13] == 'onlyintegers ':
            text += 'onlyintegers = bool(' +str(integers) +')\n'

        elif line[0:6] == 'Kspec ':
            text += 'Kspec = int(' + str(kvalue) +')\n'

        elif line[0:7] == 'Ktilde ':
            text += 'Ktilde = int(' + str(ktild) +')\n'

        elif line[0:14] == 'range_alp_def ':
            text += 'range_alp_def = int(' + str(defalp) +')\n'

        elif line[0:15] == 'points_to_calc ':
            text += 'points_to_calc = int(' + str(points) +')\n'

        elif line[0:10] == 'save_show ':
            text += 'save_show = int(' + str(svnotplot) +')\n'

        else:
            text += line

with open(config_file+'.tmp', 'w') as f:
    f.write(text)

os.rename(config_file+'.tmp',config_file)