#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 11:11:06 2017

@author: pianarol
"""

import ldf
import ldiff
import ths_dmo
import matplotlib.pyplot as plt


td = ldiff.xd
sd = ldiff.yd

#print(td)
#print(sd)
#
#print(len(td))
#print(len(sd))

#####In progress, just create the plot without analyse of which ldiff to use


###here, uses only the ldiff function
def m():
    plt.subplot(111)
    plt.loglog(ldf.t, ldf.s, c='r', label = 'drawdown')
    plt.scatter(td, sd, c='b', label = 'derivative')
    plt.grid(True)
    plt.xlabel('Time')
    plt.ylabel('Drawdown and log derivative')
    plt.title('Drawdow, Derivative, Location, Northwest')
    plt.legend()
    
    plt.show

m()

