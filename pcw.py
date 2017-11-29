#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 09:15:32 2017

@author: pianarol
"""

import hypy as hp
import scipy as sp
import numpy as np
import math
import matplotlib.pyplot as plt


###function der ###

def der(x,p):
    '''PCW_DER - Papadopulos-Cooper solution the well : Log derivative of the solution in Laplace domain
                
 Syntax: s = pcw_der( x, p)
   x[1] = Cd
   p = Laplace parameter 

 See also: pcw_lap
'''

    Sp = math.sqrt(p)
    
    k0 = sp.special.kv(0,Sp)
    
    k1 = sp.special.kv(1,Sp)
    
    s = 0.5*((2*x[0]-1)*k0**2+k1**2)/((Sp*k1+x[0]*p*k0)**2)
    
    return s

###function lap ###

def lap(x,p):
    '''PCW_LAP - Papadopulos Cooper Laplace domain solution in the well

 Syntax: s = pcw_lap( x, p)

   x[1] = Cd
   p = Laplace parameter 

 Description:
   Solution at the pumping well: 

                           K0( sqrt(p) )
  s = --------------------------------------------------------
       p [ sqrt(p)*K1( sqrt(p) ) + Cd * p * K0( sqrt(p) ) ]

 See also: pcw_dls
'''

    Sp = math.sqrt(p)
    
    k0 = sp.special.kv(0,Sp)
    
    s = k0/(p*(Sp*sp.special.kv(1,Sp)+p*x[0]*k0))

    return s

###function dls ###

def dls(x,t):
    '''PCW_DLS - Papadopulos Cooper dimensionless solution in the well

 Syntax: [s,ds] = pcw_dls( x, t)

   x(1) = Cd
   t = time 

 Description:
   Calculates the Dimensionless solution at the pumping well

    Reference: Papadopulos, I.S., and H.H.J. Cooper. 1967. Drawdown in a
    well of large diameter. Water Resources Research 3, no. 1: 241-244. 

 See also: pcw_lap
'''
    s = []
    for i in range(0, len(t)):
        s.append(hp.stefhest('pcw.lap',x,t[i]))
    
    ds = []
    for i in range(0, len(t)):
        ds.append(hp.stefhest('pcw.der',x,t[i]))
    
    
    return s,ds

###function dim ###

def dim(p,t):
    '''PCW_DIM - Papadopulos Cooper (1967) solution 

 Syntax: [s,d] = pcw_dim( p, t)

   p(1) = a  = slope of late time straight line
   p(2) = t0 = intercept of late time straight line
   p(3) = Cd = dimensionless well-bore storage coefficient

   t = time
   s = drawdonw
   d = derivative
  
% Description:
    Conputes the drawdown as a function of time with the Papadopulos and 
    Cooper (1967) solution for a pumping test in a large diameter well.
    The aquifer is confined and homogeneous. The well is fully penetrating 
    and the pumping rate constant.   
 
 
    The solution is parametrized as a function of a, to and Cd.

    The dimensionless well bore storage coefficient is:

     Cd = rc^2/(2 rw^2 S)

    a and to are the slope and time intercept of the late time straight
    line asymptotes.
     a = 0.183 Q /T
     t0 = 2.25 T t / r2 S

    NB: Note that in the original publication of Cooper et al.
    The dimensionless parameter was alpha, it is related to 
   Cd by: alpha = 1 / (2 Cd)
 
    Reference: Papadopulos, I.S., and H.H.J. Cooper. 1967. Drawdown in a
   well of large diameter. Water Resources Research 3, no. 1: 241-244. 

 Example:
   s=pcw_dim(p,t)

 See also: pcw_dmo, pcw_rpt, pcw_gss'''

    a = p[0]
    t0 = p[1]
    cd = [p[2]]
    
    var = 0.445268*t/t0
    
    
    s,d = dls(cd, var)
    
    ss = []
    dd = []
    
    for i in range(0,len(s)):
        ss.append(0.868589*a*s[i])
    
    
    for i in range(0,len(d)):
        dd.append(0.868589*a*d[i])
    
    
    
    return ss
    


###function drw ###
    
def drw():
    '''PCW_DRW - Draw the type curves of Papadopulos-Cooper (1967)

 Syntax: pcw_drw()

 Description:
   Draw a series of type curves of Papadopulos-Cooper (1967)

    Reference: Papadopulos, I.S., and H.H.J. Cooper. 1967. Drawdown in a
    well of large diameter. Water Resources Research 3, no. 1: 241-244. 

 See also: pcw_dim, pcw_dls
'''

#First figure    
    
    t = np.logspace(-2, 8)
    
    s1,ds1 = dls([10^1],t)
    s2,ds2 = dls([10^2],t)
    s3,ds3 = dls([10^3],t)
    s4,ds4 = dls([10^4],t)
    s5,ds5 = dls([10^5],t)
    
    
    st,sts = hp.ths.dls(t)
    
    fig = plt.figure()
            
    ax1 = fig.add_subplot(111)

    
    ax1.set_xlabel(r'$t_{D}$')    
    ax1.set_ylabel(r'$s_{D}$') 
    
    ax1.loglog(t,st, c='black')
    ax1.loglog(t,s1, c='b')
    ax1.loglog(t,ds1, c='b', linestyle = '--')
    ax1.loglog(t,s2, c='r')
    ax1.loglog(t,ds2, c='r', linestyle = '-.')
    ax1.loglog(t,s3, c='g')
    ax1.loglog(t,ds3, c='g', linestyle = ':')
    ax1.loglog(t,s4, c='violet')
    ax1.loglog(t,ds4, c='violet', linestyle = '--')
    ax1.loglog(t,s5, c='y')
    ax1.loglog(t,ds5, c='y', linestyle = '-.')    

    
                         
    ax1.set_ylim(ymin=1e-3)
    ax1.set_xlim(xmin=1e-2)
        
    plt.show()
    
    
#Second figure


    t = np.logspace(-1,3)
    t1 = np.power(t,10^1)
    t2 = np.power(t,10^2)
    t3 = np.power(t,10^3)
    t4 = np.power(t,10^4)
    t5 = np.power(t,10^5)
    
    s1,ds1 = dls([10^1],t1)
    s2,ds2 = dls([10^2],t2)
    s3,ds3 = dls([10^3],t3)
    s4,ds4 = dls([10^4],t4)
    s5,ds5 = dls([10^5],t5)
    
    
    
    
    fig = plt.figure()
            
    ax1 = fig.add_subplot(111)

    
    ax1.set_xlabel(r'$t_{D}$')    
    ax1.set_ylabel(r'$s_{D}$') 
    
    ax1.loglog(t,st, c='black')
    ax1.loglog(t,s1, c='b')
    ax1.loglog(t,ds1, c='b', linestyle = '--')
    ax1.loglog(t,s2, c='r')
    ax1.loglog(t,ds2, c='r', linestyle = '-.')
    ax1.loglog(t,s3, c='g')
    ax1.loglog(t,ds3, c='g', linestyle = ':')
    ax1.loglog(t,s4, c='violet')
    ax1.loglog(t,ds4, c='violet', linestyle = '--')
    ax1.loglog(t,s5, c='y')
    ax1.loglog(t,ds5, c='y', linestyle = '-.')    

    
                         
    ax1.set_ylim(ymin=1e-3)
    ax1.set_xlim(xmin=1e-2)
        
    plt.show()



###function gss ###

def gss(t,s):
    '''PCW_GSS - First guess for the parameters of the Papadopulos Cooper solution

 Syntax: p = pcw_gss(t,s)

   p(1) = a   = slope of Jacob straight line for late time
   p(2) = t0  = intercept of the Jacob straight line for late time
   p(3) = Cd  = Dimensionless coefficient (1/2alpha)

   t    = time
   s    = drawdown

 Description:
   First guess for the parameters of Papadopulos Cooper solution

 See also: pcw_dim, pcw_dmo, pcw_rpt'''

    td,d = hp.ldiffs(t,s, npoints=10)
    
    if d[-1] > 0 :
        a = math.log(10)*d[-1]
        t0 = t[-1]*math.exp(-s[-1]/d[-1])
    
    else :
        p = hp.ths.gss(t,s)
        a = p[0]
        t0 = p[1]
    
    if t0 <= 0:
        t0 = 1e-5
        
    condition = (np.greater(t,0) & np.greater(s,0))
    sp = np.extract(condition,s)
    tp = np.extract(condition,t)
    
    if not tp.all():
        print('HYTOOL: Error in pcw_gss - the vector t and s do not contain positive data')
        p = float('NaN')
        return
    else :
        cd = 0.8905356*d[-1]/sp[0]*tp[0]/t0
    
    p = []
    
    p.append(a)
    p.append(t0)
    p.append(cd)
    
    return p
    
    
###function rpt ###

def rpt(p,t,s,d, name, ttle = 'Interference test', Author = 'My name',  Rapport = 'My Rapport', filetype = 'img'):
    '''PCW_RPT - Produces the final figure and results for the Papadopulos Cooper model
    
    Syntax: pcw_rpt( p, t, s, d, ttle )
    
    p(1) = a = slope of the late time straight line 
    p(2) = t0 = intersept of late time straight line 
    p(3) = Cd = Well bore storage coefficient 
    t = measured time % s = measured drawdown 
    d(1) = Q = Pumping rate 
    d(2) = rw = Radius of well screen 
    d(3) = rc = Radius of the casing 
    ttle = Title of the figure % % Description: 
        Produces the final figure and results for the Papadopulos-Cooper model 
        
    Reference: Papadopulos, I.S., and H.H.J. Cooper. 1967. Drawdown in a 
    well of large diameter. Water Resources Research 3, no. 1: 241-244. 
    
    See also: pcw_dmo, pcw_pre, pcw_dim, pcw_gss'''

    #rename the parameters for a more intuitive check of the formulas
    q = d[0]
    rw = d[1]
    rc = d[2]
    a = p[0]
    t0 = p[1]
    cd = p[2]
    
    #Compute the transmissivity
    T = 0.1832339*q/a
    
    #Calls an internalscript that computes drawdown, derivative and residuals
    #script rpt.cmp    
    

    tc,sc,mr,sr,rms = hp.script.cmp(p,t,s,'pcw')

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
        fig.text(0.135, -0.15, 'Well radius : {:0.2g} m '.format(rw), fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.20, 'Casing radius : {:0.2g} m '.format(rc), fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.125, -0.30, 'Hydraulic parameters :', fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.35, 'Transmissivity T : {:3.1e} m²/s'.format(T), fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.125, -0.45, 'Fitting parameters :' , fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.50, 'slope a : {:0.2g} m '.format(a) , fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.55, 'intercept t0 : {:0.2g} m'.format(t0) , fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.60, 'C_D exp(2s) : {:3.1e}'.format(cd) , fontsize=14, transform=plt.gcf().transFigure)        
        fig.text(0.135, -0.66, 'mean residual : {:0.2g} m'.format(mr) , fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.70, '2 standard deviation : {:0.2g} m'.format(sr) , fontsize=14, transform=plt.gcf().transFigure)


        ax1 = fig.add_subplot(111)
        ax1.set_xlabel('Time in seconds')
        ax1.set_ylabel('Drawdown in meters')
        ax1.set_title(ttle)
        
        ax1.loglog(t, s, c='b',marker = 'o', linestyle = '', label = 'drawdown')
        ax1.loglog(td, sd, c='r',marker = 'x', linestyle = '', label = 'Derivative')
        ax1.loglog(tc, sc, c='g',   label = 'Papadopulos-Cooper (1967) Model')
        ax1.loglog(tdc, sdc, c='y',  label = 'Model derivative')
        ax1.grid(True)

        ax1.legend()
    
        plt.show()
        fig.savefig('pcw_rapport.pdf', bbox_inches = 'tight')
        
    if filetype == 'img':
        fig = plt.figure()
        fig.set_size_inches(8, 6)
        fig.text(0.125, 1, Author, fontsize=14, transform=plt.gcf().transFigure) 
        fig.text(0.125, 0.95, Rapport, fontsize=14, transform=plt.gcf().transFigure)  
        
        fig.text(0.125, -0.05, 'Test Data : ', fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.1, 'Discharge rate : {:3.2e} m³/s'.format(q), fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.15, 'Well radius : {:0.2g} m '.format(rw), fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.20, 'Casing radius : {:0.2g} m '.format(rc), fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.125, -0.30, 'Hydraulic parameters :', fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.35, 'Transmissivity T : {:3.1e} m²/s'.format(T), fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.125, -0.45, 'Fitting parameters :' , fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.50, 'slope a : {:0.2g} m '.format(a) , fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.55, 'intercept t0 : {:0.2g} m'.format(t0) , fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.60, 'C_D exp(2s) : {:3.1e}'.format(cd) , fontsize=14, transform=plt.gcf().transFigure)        
        fig.text(0.135, -0.66, 'mean residual : {:0.2g} m'.format(mr) , fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.70, '2 standard deviation : {:0.2g} m'.format(sr) , fontsize=14, transform=plt.gcf().transFigure)
       

        ax1 = fig.add_subplot(111)
        ax1.set_xlabel('Time in seconds')
        ax1.set_ylabel('Drawdown in meters')
        ax1.set_title(ttle)
        
        ax1.loglog(t, s, c='b',marker = 'o', linestyle = '', label = 'drawdown')
        ax1.loglog(td, sd, c='r',marker = 'x', linestyle = '', label = 'Derivative')
        ax1.loglog(tc, sc, c='g',   label = 'Papadopulos-Cooper (1967) Model')
        ax1.loglog(tdc, sdc, c='y',  label = 'Model derivative')
        ax1.grid(True)

        ax1.legend()
    
        plt.show()
        fig.savefig('pcw_rapport.png', bbox_inches = 'tight')        








