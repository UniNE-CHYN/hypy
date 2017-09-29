#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 25 14:57:44 2017

@author: pianarol
"""

import numpy as np
import ldf

import scipy as sp




d = 10

######################### For the time t

#make sure the array contains only float
x = np.array(ldf.t)
x = np.asfarray(x, float)
f = len(x)
a = np.log10(x[f-1])



#take only the positive values
inds_t = [i for (i, val) in enumerate(x) if val > 0]
x = [val for (i, val) in enumerate(x) if i in inds_t]

####################### For the drawdown s

#make sure the array contains only float
y = np.array(ldf.s)
y = np.asfarray(y, float)

#create w[], contains 1/standard deviation. it's use in the spline function
#give a weight to the points of the graph 
#(see documentation :https://docs.scipy.org/doc/scipy-0.19.1/reference/generated/scipy.interpolate.UnivariateSpline.html )

ysd = np.std(y)
ly = len(y)
w = []

for i in range(0, ly):
    w.append(1/ysd)






#take only the positive values
inds_s = [i for (i, val) in enumerate(y) if val > 0]
y = [val for (i, val) in enumerate(y) if i in inds_s]

f = len(x)

#xi and yi
#changing k,s and w affects the graph, make it better or worse

xi = np.logspace(np.log10(x[0]), np.log10(x[f-1]),  num = f, endpoint = True, base = 10.0, dtype = np.float64)

spl = sp.interpolate.UnivariateSpline(x,y,w, k = 5, s = 0.0099)
yi = spl(xi)


#xd and yd

xd = xi[1:len(xi)-1]
yd = xd*(yi[2:len(yi)+1]-yi[0:len(yi)-2])/(xi[2:len(xi)+1]-xi[0:len(xi)-2])



