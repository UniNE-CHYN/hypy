#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 09:42:40 2017

@author: pianarol
"""

import ldf
import numpy as np
import math


d = 2

x = ldf.t
y = ldf.s

######################### For the time t

#make sure the array contains only float
x = np.array(ldf.t)
x = np.asfarray(x, float)


#take only the positive values
inds_t = [i for (i, val) in enumerate(x) if val > 0]
x = [val for (i, val) in enumerate(x) if i in inds_t]


###transform the value of the array X into logarithms 
logx = []

for i in range(0, len(x)):
    logx.append(math.log(x[i]))


dx = np.diff(logx)
dx1 = dx[0:len(dx)-2*d+1]
dx2 = dx[3:len(dx)] ###########attention au 3, normalement c'est 1*d, mais ne donne pas le résultat attendu






####################### For the drawdown s

#make sure the array contains only float
y = np.array(ldf.s)
y = np.asfarray(y, float)

#take only the positive values
inds_s = [i for (i, val) in enumerate(y) if val > 0]
y = [val for (i, val) in enumerate(y) if i in inds_s]


dy = np.diff(y)
dy1 = dy[0:len(dx)-2*d+1]
dy2 = dy[3:len(dy)]


#xd and yd

xd = x[2:len(x)-2]
yd = (dx2*dy1/dx1+dx1*dy2/dx2)/(dx1+dx2)
