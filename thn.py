#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 23 16:07:07 2017

@author: pianarol
"""

import hypy as hp
import scipy as sp
import numpy as np
import math
import matplotlib.pyplot as plt

###function dim ###

def dim(p,t):
    '''THN_DIM - Theis model with a no-flow boundary

 Syntax: s = thn_dim(p,t)
   p(1) = a  = slope of Jacob Straight Line
   p(2) = t0 = intercept of the first segment of straight line
   p(3) = ti = time of intersection between the 2 straight lines

 Description:
   Calculate the drawdown at time t for confined aquifer with a no-flow
   boundary

 Example:
   s = thn_dim( p,t )

 See also: thn_dmo, thn_dim, thn_gss'''

    a = p[0]
    
    s = hp.ths.dim([a,p[1]],t) + hp.ths.dim([a,p[2]],t)
    
    return s


###function dls ###
    
def dls(p,t):
    '''THN_DLS - Dimensionless drawdown of the Theis model with a no-flow boundary

 Syntax: s = thn_dls(p,t)

   p(1) = r1/r2  
   with r1 = radius to the pumping well 
   and  r2 = radius to the image well 

   provides the dimensionless drawdown at reduced time t'''

    t = 0.25/t
    
    s = 0.5*sp.special.expn(1,t)+0.5*sp.special.expn(1,t*p[0]**2)
    
    return s    


###function drw ###

def drw():
    '''THN_DRW - Type curves of the Theis model with a no-flow boundary

 Syntax: thn_drw()

 Description:
   Draw a series of type curves of Theis (1935) model with a no-flow 
   boundary

 See also: thn_dim, thn_dls'''
 
#First figure    
    
    td = np.logspace(-1, 4)
    
    R1 = [1.3,3.3,10,30]
    R2 = [3.3,1.3]
    R3 = [10,0]
    R4 = [30,0]
    
    s1 = hp.thn.dls(R1,td)
    s2 = hp.thn.dls(R2,td)
    s3 = hp.thn.dls(R3,td)
    s4 = hp.thn.dls(R4,td)
    xd1,yd1 = hp.ldiff(td,s1)
    xd2,yd2 = hp.ldiff(td,s2)
    xd3,yd3 = hp.ldiff(td,s3)
    xd4,yd4 = hp.ldiff(td,s4)    
    
    tds,tdss = hp.ths.dls(td)
    
    
    fig = plt.figure()
            
    ax1 = fig.add_subplot(111)

    
    ax1.set_xlabel(r'$t_{D}/r²_{D}$')    
    ax1.set_ylabel(r'$s_{D}$') 
    
    ax1.loglog(td,s1, c='b')
    ax1.loglog(xd1,yd1, c='b', linestyle = '--')
    ax1.loglog(td,s2, c='r')
    ax1.loglog(xd2,yd2, c='r', linestyle = '-.')
    ax1.loglog(td,s3, c='g')
    ax1.loglog(xd3,yd3, c='g', linestyle = ':')
    ax1.loglog(td,s4, c='y')
    ax1.loglog(xd4,yd4, c='y', linestyle = '--')
    
    ax1.loglog(td,tds, c = 'black')
    
#    ax1.legend()                          
    ax1.set_ylim(ymin=1e-2)
    ax1.set_xlim(xmin=1e-1)
        
    plt.show()
    
    
#Second figure

    Rd = [33]

    s = hp.thn.dls(Rd,td)

    xd,yd = hp.ldiff(td,s)
    tdj = hp.jcb.dls(td)
    
    tdrd = []
    
    for i in range(0,len(tdj)):
        tdrd.append(2*tdj[i])
    
    fig2 = plt.figure()
            
    ax2 = fig2.add_subplot(111)

    ax2.semilogx(td,s,marker = 'd', linestyle = '', c = 'b')
    ax2.semilogx(xd,yd, marker = 'd', linestyle ='', c = 'r')
    ax2.semilogx(td,tdj, c='g', linestyle = '--')
    ax2.semilogx(td,tdrd,c = 'black', linestyle = '--')
    
    ax2.set_xlabel('t')
    ax2.set_ylabel('s')
    
    ax2.set_ylim(ymin = 0)

    plt.show()             


###function gss ###

def gss(t,s):
    '''THN_GSS - First guess for the parameters of the Theis model with a no-flow boundary

 Syntax: p = thn_gss(t,s)

   p(1) = a  = slope of Jacob straight line 
   p(2) = t0 = intercept of the Jacob straight line 
   p(3) = ti = time of intersection between the 2 straight lines
   t    = time
   s    = drawdown

 Description:
   First guess for the parameters of theis solution with a no-flow
   boundary

 See also: thn_dmo, thn_rpt, thn_dim
'''

    #Automatic identification of the "control" points
    td,d = hp.ldiffs(t,s, npoints=10) #First log derivative
    tdd,dd = hp.ldiffs(td,d, npoints=10) #Second log derivative

    
    #Calculation of the parameters of the model
    tmp = np.amax(dd)

    tmp = np.float(tmp)
    i = np.argmax(dd)

    ti = tdd[i-1]
    
    #Slope of Jacob's straight line
    
    a = d[-1]*2.30/2
    
    #Origin of jacob straight line
    t0 = math.pow(10,(a*math.log10(t[-1]*t[-1]/ti)-s[-1])/a)
    
    p = []
    
    p.append(a)
    p.append(t0)
    p.append(ti)
    
    return p    


###function rpt ###
    
def rpt(p,t,s,d, name, ttle = 'Interference test', Author = 'My name',  Rapport = 'My Rapport', filetype = 'img'):
    '''THN_RPT - Produces the final figure and results for the Theis model with a no flow boundary
    
    Syntax: thn_rpt( p, t, s, d, ttle )
    
    p = parameters of the model 
    t = measured time 
    s = measured drawdown 
    d(1) = Q = Pumping rate 
    d(2) = r = Distance between the observation and the pumping well 
    ttle = Title of the figure (Optional) 
    
    Description: 
         Produces the final figure and results for Theis model (1935) with 
          a no flow boundary. 
    See also: thn_dmo, thn_dim, thn_gss'''
        
    #rename the parameters for a more intuitive check of the formulas
    a = p[0]
    t0 = p[1]
    ti = p[2]
    q = d[0]
    r = d[1]
    
    #Compute the transmissivity, storativity and radius of influence
    T = 0.1832339*q/a
    S = 2.458394*T*t0/r**2
    Ri = math.sqrt(2.2458394*T*ti/S)    

    #Calls an internalscript that computes drawdown, derivative and residuals
    #script rpt.cmp    
    

    tc,sc,mr,sr,rms = hp.script.cmp(p,t,s,'thn')
    
    #script rpt_plt
    
    #calculate the derivative of the data
    td, sd = hp.ldiffs(t,s, npoints=30)
    #keep only positive derivatives
    td, sd = hp.hyclean(td,sd)
    
    #compute the derivative of the model
    tdc,sdc = hp.ldiff(tc,sc)
    #keep only positive derivatives
    tdc,sdc = hp.hyclean(tdc,sdc)
    
    #plots the data and model in bi-logarithmic scale
    if filetype == 'pdf':
        
        fig = plt.figure()
        fig.set_size_inches(8, 6)
        fig.text(0.125, 1, Author, fontsize=14, transform=plt.gcf().transFigure) 
        fig.text(0.125, 0.95, Rapport, fontsize=14, transform=plt.gcf().transFigure)

        fig.text(0.125, -0.05, 'Test Data : ', fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.1, 'Discharge rate : {:3.2e} m³/s'.format(q), fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.15, 'Radial distance : {:0.2g} m '.format(r), fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.125, -0.25, 'Hydraulic parameters :', fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.3, 'Transmissivity T : {:3.1e} m²/s'.format(T), fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.35, 'Storativity S : {:3.1e} '.format(S), fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.40, 'Distance to image well Ri : {:0.2g} m'.format(Ri) , fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.125, -0.5, 'Fitting parameters :' , fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.55, 'slope a : {:0.2g} m '.format(a) , fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.60, 'intercept t0 : {:0.2g} m'.format(t0) , fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.65, 'intercept ti : {:0.2g} m'.format(ti) , fontsize=14, transform=plt.gcf().transFigure)        
        fig.text(0.135, -0.70, 'mean residual : {:0.2g} m'.format(mr) , fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.75, '2 standard deviation : {:0.2g} m'.format(sr) , fontsize=14, transform=plt.gcf().transFigure)


        ax1 = fig.add_subplot(111)
        ax1.set_xlabel('Time in seconds')
        ax1.set_ylabel('Drawdown in meters')
        ax1.set_title(ttle)
        
        ax1.loglog(t, s, c='r',marker = '+', linestyle = '', label = 'drawdown')
        ax1.loglog(td, sd, c='b',marker = 'x', linestyle = '', label = 'Derivative')
        ax1.loglog(tc, sc, c='g',   label = 'Theis (1935) Model')
        ax1.loglog(tdc, sdc, c='y',  label = 'Model derivative')
        ax1.grid(True)

        ax1.legend()
    
        plt.show()
        fig.savefig('thn_rapport.pdf', bbox_inches = 'tight')
        
    if filetype == 'img':
        fig = plt.figure()
        fig.set_size_inches(8, 6)
        fig.text(0.125, 1, Author, fontsize=14, transform=plt.gcf().transFigure) 
        fig.text(0.125, 0.95, Rapport, fontsize=14, transform=plt.gcf().transFigure)  
        
        fig.text(0.125, -0.05, 'Test Data : ', fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.1, 'Discharge rate : {:3.2e} m³/s'.format(q), fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.15, 'Radial distance : {:0.2g} m '.format(r), fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.125, -0.25, 'Hydraulic parameters :', fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.3, 'Transmissivity T : {:3.1e} m²/s'.format(T), fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.35, 'Storativity S : {:3.1e} '.format(S), fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.40, 'Distance to image well Ri : {:0.2g} m'.format(Ri) , fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.125, -0.5, 'Fitting parameters :' , fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.55, 'slope a : {:0.2g} m '.format(a) , fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.60, 'intercept t0 : {:0.2g} m'.format(t0) , fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.65, 'intercept ti : {:0.2g} m'.format(ti) , fontsize=14, transform=plt.gcf().transFigure)        
        fig.text(0.135, -0.70, 'mean residual : {:0.2g} m'.format(mr) , fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.75, '2 standard deviation : {:0.2g} m'.format(sr) , fontsize=14, transform=plt.gcf().transFigure)        

        ax1 = fig.add_subplot(111)
        ax1.set_xlabel('Time in seconds')
        ax1.set_ylabel('Drawdown in meters')
        ax1.set_title(ttle)
        
        ax1.loglog(t, s, c='r',marker = '+', linestyle = '', label = 'drawdown')
        ax1.loglog(td, sd, c='b',marker = 'x', linestyle = '', label = 'Derivative')
        ax1.loglog(tc, sc, c='g',   label = 'Theis (1935) Model')
        ax1.loglog(tdc, sdc, c='y',  label = 'Model derivative')
        ax1.grid(True)

        ax1.legend()
    
        plt.show()
        fig.savefig('thn_rapport.png', bbox_inches = 'tight')    
