#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 15:33:01 2017

@author: pianarol
"""

import numpy as np
import ldf
import math

#LDIFF - Approximate logarithmic derivative with centered differences

######################### For the time t

#make sure the array contains only float
x = np.array(ldf.t)
x = np.asfarray(x, float)


#take only the positive values
inds_t = [i for (i, val) in enumerate(x) if val > 0]
x = [val for (i, val) in enumerate(x) if i in inds_t]

#Calculate the difference
dx = np.diff(x)



####################### For the drawdown s

#make sure the array contains only float
y = np.array(ldf.s)
y = np.asfarray(y, float)

#take only the positive values
inds_s = [i for (i, val) in enumerate(y) if val > 0]
y = [val for (i, val) in enumerate(y) if i in inds_s]

#Calculate the difference
dy = np.diff(y)

####calculate the point

#xd 
xd = []

for i in range(0, len(x)-1):
    xd.append( math.sqrt(x[i]*x[i+1]))


#yd
yd = xd*(dy/dx)








    
