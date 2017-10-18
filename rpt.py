#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 14:07:08 2017

@author: pianarol
"""

import numpy as np
import math
import hypy as hp



###function cmp ###
# RPT_CMP
#
# Description: 
#   Internal script used for the ***_rpt functions.
#   It computes the drawdown, derivative and residuals of a given model.
#
# See also: rpt_lgd, rpt_llt, rpt_plt


def cmp(p,t,s,name):
    #keep only the positive time
    t,s = hp.hyclean(t,s)
    #define regular points to plot the calculated drawdown
    tc = np.logspace(np.log10(t[0]), np.log10(t[len(t)-1]),  num = len(t), endpoint = True, base = 10.0, dtype = np.float64)
    
    #compute the drawdown with the model
    if name == 'ths':
        sc = hp.ths.dim(p,tc)
    if name == 'Del' : 
        sc = hp.Del.dim(p,tc)
    #keep only the positive drawdown
    tc,sc = hp.hyclean(tc,sc)
    
    #Compute the residuals and their statistics
    residuals = s - hp.ths.dim(p,t)
    mr = np.mean(residuals)
    sr = 2 * np.nanstd(residuals)
    rms = math.sqrt(np.mean(residuals**2))

    return tc,sc,mr,sr,rms


