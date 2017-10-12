#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  9 13:59:33 2017

@author: pianarol
"""

import numpy as np
import math
import hypy as hp

###function gss ###

#JCB_GSS  - Linear Least Square fitting of Jacob solution
#
# Syntax: p = jcb_gss(t,s,tmin,tmax)
#
#   p(1) = a  = slope of Jacob straight line
#   p(2) = t0 = intercept of the Jacob straight line
#
#   t    = time 
#   s    = drawdown
#   tmin = optional argument: start of the fitting period
#   tmax = optional argument: end   of the fitting period
#
# Description:
#   First guess for the parameters of Jacob solution
#
# See also: jcb_dim, jcb_dmo, jcb_rpt

def gss(t,s):
    tmin = t[0]
    tmax = t[len(t)-1]
    
    
    
    
    logt = []
    for i in range(0, len(t)):
        logt.append(math.log10(t[i]))
    t1 = []
    for i in range(0, len(logt)):
        t1.append(1)
    
    g = np.transpose([logt, t1])
    g1 = np.transpose(g)
    p1 = np.linalg.inv(np.dot(g1,g))
    p2 = np.dot(p1,g1)

    p = np.dot(p2,s)
    
    a = p[0]
    c = p[1]
    t0 = math.pow(10, -c/a)
    p[1] = t0
    return a, t0
    
###function dls ###
#JCB_DLS - Jacob dimensionless drawdown
#
# Syntax: [s,d] = jcb_dls(td)
#
# Description:
#   provides the dimensionless drawdown and derivative at reduced time td
#
# See also: jcb_dmo


def dls(td):
    s = []
    
    for i in range(0, len(td)):
        s.append(0.5*(math.log(4*td[i])-0.5772))
    return s
    


