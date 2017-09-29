#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 11:11:06 2017

@author: pianarol
"""

##########Still in progress, need to contain a loop to chose which function need to be used

import ldf
import ldiff
import ldiffs
import ldiffb
import ldiffh

import matplotlib.pyplot as plt


###here, uses only the ldiff function
def ldiff_plot():
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Drawdown and log derivative')
    ax1.set_title('Drawdow, Derivative, Location, Northwest')
    ax1.loglog(ldf.t, ldf.s, c='r', label = 'drawdown')
    ax1.scatter(ldiff.xd, ldiff.yd, c='b', label = 'derivative')
    ax1.grid(True)

    ax1.legend()
    
    plt.show

if __name__ == "__main__":              #make the function works, no idea why
    ldiff_plot()


####here, uses only the ldiffs function
def ldiffs_plot():
    fig = plt.figure()
    ax2 = fig.add_subplot(111)
    ax2.set_xlabel('Time')
    ax2.set_ylabel('Drawdown and log derivative')
    ax2.set_title('Drawdow, Derivative, Location, Northwest')
    ax2.loglog(ldf.t, ldf.s, c='r', label = 'drawdown')
    ax2.scatter(ldiffs.xd, ldiffs.yd, c='b', label = 'derivative')
    ax2.grid(True)

    ax2.legend()
    
    plt.show

if __name__ == "__main__":              #make the function works, no idea why
    ldiffs_plot()

####here, uses only the ldiffb function
def ldiffb_plot():
    fig = plt.figure()
    ax2 = fig.add_subplot(111)
    ax2.set_xlabel('Time')
    ax2.set_ylabel('Drawdown and log derivative')
    ax2.set_title('Drawdow, Derivative, Location, Northwest')
    ax2.loglog(ldf.t, ldf.s, c='r', label = 'drawdown')
    ax2.scatter(ldiffb.xd, ldiffb.yd, c='b', label = 'derivative')
    ax2.grid(True)

    ax2.legend()
    
    plt.show

if __name__ == "__main__":              #make the function works, no idea why
    ldiffb_plot()

####here, uses only the ldiffh function
def ldiffh_plot():
    fig = plt.figure()
    ax2 = fig.add_subplot(111)
    ax2.set_xlabel('Time')
    ax2.set_ylabel('Drawdown and log derivative')
    ax2.set_title('Drawdow, Derivative, Location, Northwest')
    ax2.loglog(ldf.t, ldf.s, c='r', label = 'drawdown')
    ax2.scatter(ldiffh.xd, ldiffh.yd, c='b', label = 'derivative')
    ax2.grid(True)

    ax2.legend()
    
    plt.show

if __name__ == "__main__":              #make the function works, no idea why
    ldiffh_plot()


