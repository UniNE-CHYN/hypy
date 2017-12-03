# -*- coding: utf-8 -*-
"""
Created on Sun Dec  3 14:10:24 2017

@author: Loïc
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 27 14:19:25 2017

@author: pianarol
"""

import hypy as hp
import scipy as sp
import numpy as np
import math
import matplotlib.pyplot as plt


   

###function lap ###

def lap(x,p):
    '''BLT_LAP - Boulton (1963) Laplace dimensionless domain solution 

 Syntax: s = hp.blt.lap( x, p)

   x[1] = sigma = S / S_y
   x[2] = phi   = ( alpha r^2 S ) / T 
   p    = Laplace parameter 

 See also: blt_dls'''
 
    s = sp.special.kv(0, math.sqrt(p+x[1]*p/(x[0]*(p+x[1]))))/p
    
    return s
    
    
###function dls ###

def dls(x,t):
    '''BLT_DLS - Dimensionless drawdown of the Boulton model for an unconfined aquifer. 

 Syntax: s = hp.blt.dls( x, t)

   Dimensionless solution 
        
   x[1] = sigma = S / S_y
   x[2] = phi   = ( alpha r^2 S ) / T 
   t    = dimensionless time
 See also: blt_lap''' 

    s = []
    for i in range(0, len(t)):
        s.append(hp.stefhest('blt.lap',x,t[i])) 
        
    return s
    
    
###function dim ###

def dim(p,t):
    '''BLT_DIM(p,t) Boulton (1963) model for unconfined aquifer

   Syntax: s = hp.blt.dim(p,t)

       p = vector of the parameters
        
       p(1) = a : slope of the late time Jacob's straight line in m
       p(2) = t0 : intercept with the horizontal axis for s = 0
       p(3) = t1
       p(4) = phi
       t    = dimensionless time 

       t = time in seconds
       s = drawdown in meters


   Description:
    Computes drawdon in the aquifer for a constant rate pumping test.
    The aquifer is unconfined and homogeneous.

   See also: blt_rpt, blt_gss'''

    td = 0.445268*t/p[1]
    
    sds = dls([p[1]/(p[1]+p[2]),2*p[3]*p[1]],td)
    
    sd = []
    
    for i in range(0,len(sds)):
        sd.append(0.868589*p[0]*sds[i])

    return sd    
    
    
###function drw ###

###Not sure if useful and coded right

#def drw():
#    '''BLT_DRW - Draw the type curves of Boulton (1963)
#
# Syntax: hp.blt.drw()
#
# Description:
#   Draw a series of type curves of Boulton (1963).
#
# See also: blt_dim, blt_dls'''    
#    
##First figure    
#    
#    sigma = 0.01
#    t = np.logspace(-2, 6)
#    
#    s1 = dls([sigma,10^1],t)
#    s2 = dls([sigma,10^2],t)
#    s3 = dls([sigma,10^3],t)
#    s4 = dls([sigma,10^4],t)
#    
#
#    
#    st1,sts = hp.ths.dls(t)
#    st2 = 0.5*sp.special.expn(1,1+sigma)/(4*t*sigma)
#    
#    fig = plt.figure()
#            
#    ax1 = fig.add_subplot(111)
#
#    
#    ax1.set_xlabel(r'$t_{D}$')    
#    ax1.set_ylabel(r'$s_{D}$') 
#    
#    ax1.loglog(t,st1, c='black')
#    ax1.loglog(t,st2,c='orange')
#    ax1.loglog(t,s1, c='blue')
#
#    ax1.loglog(t,s2, c='red')
#
#    ax1.loglog(t,s3, c='green')
# 
#    ax1.loglog(t,s4, c='violet')
#  
#
#    
#                         
#    ax1.set_ylim(ymin=1e-3)
#    ax1.set_xlim(xmin=1e-2)
#        
#    plt.show()        
    



###function gss ###

def gss(t,s):
    '''BLT_GSS - First guess for the parameters of the Boulton model

 Syntax:  p = hp.blt.gss(t,s)
   p(1) = a  = slope of Jacob straight line for late time
   p(2) = t0 = intercept with the horizontal axis for s = 0 for the early
               asymptote
   p(3) = t1 = intercept with the horizontal axis for s = 0 for the late
               asymptote
   p(4) = phi

   t    = time
   s    = drawdown

 Description:
   First guess for the parameters of Boulton solution.

 See also: blt_dmo, blt_rpt, blt_dim'''  

    p = [1,2,3,4]   
    
    p[1] = t[0]
    
    n = round(len(t)/3)
    
    
    t1 = t[n-1:len(t)]
    s1 = s[n-1:len(s)]
    
    
    
    pj = hp.jcb.gss(t1,s1)
    

    p[3] = 1e-4
    p[0] = pj[0]
    p[2] = pj[1]
    
    
    return p    
    
    
    
###function rpt ###

def rpt(p,t,s,d, name, ttle = 'Interference test', Author = 'My name',  Rapport = 'My Rapport', filetype = 'img'):
    '''BLT_RPT - Reports the results of the interpretation with the Boulton model
    
    Syntax: hp.blt.rpt( p, t, s, d, ttle )
    
    p = parameters of the model 
    p(1) = a = slope of the straight lines 
    p(2) = t0 = intercept with the horizontal axis the first line 
    p(3) = t1 = intercept with the horizontal axis the second line 
    p(4) = phi = empirical parameter that trigger the delay 
    d(1) = q = pumping rate 
    d(2) = r = distance between pumping well and piezometer 
    t = measured time 
    s = measured drawdown 
    ttle = Title of the figure 
    
    Description: 
        Produces the final figure and results for the Boulton model (1963). 
        See also: blt_dmo, blt_dim, blt_pre, blt_gss'''
        
    #Rename the parameters for a more intuitive check of the formulas
    a = p[0]
    t0 = p[1]
    t1 = p[2]
    phi = p[3]
    
    q = d[0]
    r = d[1]

    #Compute the transmissivity, storativity and radius of influence
    T=0.1832339*q/a
    S=2.2458394*T*t0/r**2 
    omegad=2.2458394*T*t1/r**2-S
    Ri=2*math.sqrt(T*t[-1]/omegad)

    #Calls an internalscript that computes drawdown, derivative and residuals
    #script rpt.cmp    
    

    tc,sc,mr,sr,rms = hp.script.cmp(p,t,s,'blt')

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
        fig.text(0.135, -0.30, 'Transmissivity T : {:3.1e} m²/s'.format(T), fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.35, 'Storativity S : {:3.1e}'.format(S), fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.40, 'Drainage Porosity : {:3.1e}'.format(omegad), fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.45, 'Radius of Investigation Ri : {:0.2g} m'.format(Ri), fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.125, -0.55, 'Fitting parameters :' , fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.60, 'slope a : {:0.2g} m '.format(a) , fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.65, 'intercept t0 : {:0.2g} s'.format(t0) , fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.70, 'intercept t1 : {:0.2g} s'.format(t1) , fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.75, 'phi : {:0.2g} m'.format(phi) , fontsize=14, transform=plt.gcf().transFigure)        
        fig.text(0.135, -0.80, 'mean residual : {:0.2g} m'.format(mr) , fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.85, '2 standard deviation : {:0.2g} m'.format(sr) , fontsize=14, transform=plt.gcf().transFigure)


        ax1 = fig.add_subplot(111)
        ax1.set_xlabel('Time in seconds')
        ax1.set_ylabel('Drawdown in meters')
        ax1.set_title(ttle)
        
        ax1.loglog(t, s, c='b',marker = 'o', linestyle = '', label = 'drawdown')
        ax1.loglog(td, sd, c='r',marker = 'x', linestyle = '', label = 'Derivative')
        ax1.loglog(tc, sc, c='g',   label = 'Boulton (1963) Model')
        ax1.loglog(tdc, sdc, c='y',  label = 'Model derivative')
        ax1.grid(True)

        ax1.legend()
    
        plt.show()
        fig.savefig('blt_rapport.pdf', bbox_inches = 'tight')
        
    if filetype == 'img':
        fig = plt.figure()
        fig.set_size_inches(8, 6)
        fig.text(0.125, 1, Author, fontsize=14, transform=plt.gcf().transFigure) 
        fig.text(0.125, 0.95, Rapport, fontsize=14, transform=plt.gcf().transFigure)  
        
        fig.text(0.125, -0.05, 'Test Data : ', fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.1, 'Discharge rate : {:3.2e} m³/s'.format(q), fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.15, 'Radial distance : {:0.2g} m '.format(r), fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.125, -0.25, 'Hydraulic parameters :', fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.30, 'Transmissivity T : {:3.1e} m²/s'.format(T), fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.35, 'Storativity S : {:3.1e}'.format(S), fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.40, 'Drainage Porosity : {:3.1e}'.format(omegad), fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.45, 'Radius of Investigation Ri : {:0.2g} m'.format(Ri), fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.125, -0.55, 'Fitting parameters :' , fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.60, 'slope a : {:0.2g} m '.format(a) , fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.65, 'intercept t0 : {:0.2g} s'.format(t0) , fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.70, 'intercept t1 : {:0.2g} s'.format(t1) , fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.75, 'phi : {:0.2g} m'.format(phi) , fontsize=14, transform=plt.gcf().transFigure)        
        fig.text(0.135, -0.80, 'mean residual : {:0.2g} m'.format(mr) , fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.85, '2 standard deviation : {:0.2g} m'.format(sr) , fontsize=14, transform=plt.gcf().transFigure)
       

        ax1 = fig.add_subplot(111)
        ax1.set_xlabel('Time in seconds')
        ax1.set_ylabel('Drawdown in meters')
        ax1.set_title(ttle)
        
        ax1.loglog(t, s, c='b',marker = 'o', linestyle = '', label = 'drawdown')
        ax1.loglog(td, sd, c='r',marker = 'x', linestyle = '', label = 'Derivative')
        ax1.loglog(tc, sc, c='g',   label = 'Boulton (1963) Model')
        ax1.loglog(tdc, sdc, c='y',  label = 'Model derivative')
        ax1.grid(True)

        ax1.legend()
    
        plt.show()
        fig.savefig('blt_rapport.png', bbox_inches = 'tight')    
         