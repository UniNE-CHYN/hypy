#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  5 15:55:54 2017

@author: pianarol
"""

import scipy as sp
import math
import numpy as np
import hypy as hp
import matplotlib.pyplot as plt


###Contains all the functions of the module ths####


####function dim ###



def dim(p,t):
    '''THS_DIM - Compute drawdown with the Theis (1935) solution

 Syntax: s = ths_dim( p, t) 

   p(1) = a  = slope of Jacob Straight Line in meters
   p(2) = t0 = intercept with the horizontal axis for s = 0
   t = measured time
   s = measured drawdown

 Description:
   The Theis (1935) solution assumes that the aquifer is confined,
   homogeneous, isotropic, and infinite. The well as radius that is
   negligible. It is pumped at a constant rate Q and is 100 percent
   efficient.   

   Under these assumptions, Theis solution can be expressed as: 

       s(r,t) = Q/(4 pi T) E1( r2S / 4Tt)

   where Q is the pumping rate, T the transmissivity, r the radial 
   distance between the pumping well and the observation well, 
   S the storativity coefficient and t the time.  

   To interpret field data with the Theis solution, it is expressed as a
   function of two fitting parameters a and to which are defined as:   

   a = 0.183 Q /T
   t0 =  r2 S / 2.25 T 

 Example:
   s = ths_dim( p,t )

 See also: ths_dmo, ths_gss, ths_rpt'''
    td = []
    for i in range(0, len(t)):
        td.append(0.5628*p[1]/t[i])
        
    s = p[0]/math.log(10)*sp.special.expn(1,td)
    
    tdn = [ -x for x in td]
    d = p[0]/math.log(10)*np.exp(tdn)

    return s

####function dls ###





def dls(td):
    '''THS_DLS - Dimensionless drawdown of the Theis model

 Syntax: [sd,dd] = ths_dls(td)

 Description:
   Calculates the dimensionless drawdown sd and the dimensionless 
   derivative dd for a given dimensionless reduced time td/rd^2

 See also: ths_lap'''
    td2 = []
    dd = []
    for i in range(0, len(td)):
        td2.append((0.25/td[i]))
        dd.append(0.5*np.exp(-1.0*td2[i]))

    sd = 0.5*sp.special.expn(1,td2)
    
    
    return sd, dd

####function drw ###


def drw():
    '''THS_DRW - Type curves of the Theis model.

 Syntax: ths_drw()

 Description:
   Draw a series of type curves of Theis (1935)

 See also: ths_dim, ths_dls
'''

    td = np.logspace(-1, 4)
    sd, dd = dls(td)
    fig = plt.figure()
            
    ax1 = fig.add_subplot(111)

       
    ax1.set_xlabel(r'$t_{D}/r²_{D} = T t / S r²$')
    ax1.set_ylabel(r'$s_{D} = 2*\pi*s/Q$')    
                                                 
    ax1.loglog(td,sd, c='b', label='Theis') 
    ax1.plot(td,dd, c='r', linestyle = '--', label='Derivative')        
    ax1.legend()                          
            
    plt.show()                              

    fig2 = plt.figure()
    sj = hp.jcb.dls(td)        
    ax2 = fig2.add_subplot(111)

        
    ax2.set_xlabel(r'$t_{D}/r²_{D}$')    
    ax2.set_ylabel(r'$s_{D}$')     
    
    ax2.semilogx(td,sd, c='b', linestyle = '--', label='Theis') 
    ax2.plot(td,sj, c='r', label='Jacob')        
    ax2.legend()                           
            
    plt.show()                              

### function gss ###



def gss(t,s):
    '''THS_GSS - First guess for the parameters of the Theis model.

 Syntax:  p = ths_gss(t,s)

   p(1) = a  = slope of Jacob straight line for late time
   p(2) = t0 = intercept with the horizontal axis for s = 0
   t    = time
   s    = drawdown

 Description:
   First guess for the parameters of theis solution

 See also: ths_dim, ths_dmo, ths_rpt, ezwt'''

    if np.shape(t) == 1:
        t = np.transpose(t)
        s = np.transpose(s) #contrôler si c'est ce qui est souhaité
    
    n = round(len(t)/3)
    t = t[n:len(t)]
    s = s[n:len(s)]
    p = hp.jcb.gss(t,s)
    
    return p






### function jac ###


def jac(p,t):
    '''THS_JAC - Jacobian matrix of the Theis function

 Syntax: j = ths_jac( p, t)

    j(1,:) = ds / dp(1) 
    j(2,:) = ds / dp(2)
'''
    td = []
    for i in range(0, len(t)):
        td.append(0.5625*p[1]/t[i])
    j1 = []
    for i in range(0, len(td)):
        j1.append(sp.special.expn(1,td[i])/math.log(10))
       
    tdn = [-x for x in td]     
    j2 = []
    for i in range(0, len(td)):
        j2.append(p[0]*np.exp(tdn[i])/math.log(10)/p[1])

    j = [j1,j2]   

    return j 

###unfction rpt ###



def rpt(p,t,s,d, name, ttle = 'Interference test', Author = 'My name',  Rapport = 'My Rapport', filetype = 'img'):
    '''THS_RPT - Reports graphically the results of a pumping test interpreted with the Theis (1935) model. 
 Syntax: ths_rpt( p, t, s, d, name, ttle, Author, Rapport, filetype )
 p = parameters of the model 
 p(1) = a = slope of the straight line
 p(2) = t0 = intercept with the horizontal axis for s = 0                                                                                        
 t = measured time % s = measured drawdown 
 d(1) = q = Pumping rate 
 d(2) = r = distance from the pumping well 
 ttle = Title of the figure 
 Description: 
     Produces the final figure and results for Theis model (1935).
See also: ths_dmo, ths_dim, ths_gss'''
    
    #rename the parameters for a more intuitive check of the formulas
    a = p[0]
    t0 = p[1]
    q = d[0]
    r = d[1]
    
    #Compute the transmissivity, storativity and radius of influence
    T = 0.1832339*q/a
    S = 2.458394*T*t0/r**2
    Ri = 2*math.sqrt(T*t[len(t)-1]/S)

    
    #Calls an internalscript that computes drawdown, derivative and residuals
    #script rpt.cmp    
    
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
        fig.text(0.135, -0.3, 'Transmissivity T : {:3.2e} m²/s'.format(T), fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.35, 'Storativity S : {:3.2e} '.format(S), fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.40, 'Radius of Investigation Ri : {:0.2g} m'.format(Ri) , fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.125, -0.5, 'Fitting parameters :' , fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.55, 'slope a : {:0.2g} m '.format(a) , fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.60, 'intercept t0 : {:0.2g} m'.format(t0) , fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.65, 'mean residual : {:0.2g} m'.format(mr) , fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.70, '2 standard deviation : {:0.2g} m'.format(sr) , fontsize=14, transform=plt.gcf().transFigure)


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
        fig.savefig('ths_rapport.pdf', bbox_inches = 'tight')
        
    if filetype == 'img':
        fig = plt.figure()
        fig.set_size_inches(8, 6)
        fig.text(0.125, 1, Author, fontsize=14, transform=plt.gcf().transFigure) 
        fig.text(0.125, 0.95, Rapport, fontsize=14, transform=plt.gcf().transFigure)  
        
        fig.text(1, 0.85, 'Test Data : ', fontsize=14, transform=plt.gcf().transFigure)
        fig.text(1.05, 0.8, 'Discharge rate : {:3.2e} m³/s'.format(q), fontsize=14, transform=plt.gcf().transFigure)
        fig.text(1.05, 0.75, 'Radial distance : {:0.2g} m '.format(r), fontsize=14, transform=plt.gcf().transFigure)
        fig.text(1, 0.65, 'Hydraulic parameters :', fontsize=14, transform=plt.gcf().transFigure)
        fig.text(1.05, 0.6, 'Transmissivity T : {:3.2e} m²/s'.format(T), fontsize=14, transform=plt.gcf().transFigure)
        fig.text(1.05, 0.55, 'Storativity S : {:3.2e} '.format(S), fontsize=14, transform=plt.gcf().transFigure)
        fig.text(1.05, 0.5, 'Radius of Investigation Ri : {:0.2g} m'.format(Ri) , fontsize=14, transform=plt.gcf().transFigure)
        fig.text(1, 0.4, 'Fitting parameters :' , fontsize=14, transform=plt.gcf().transFigure)
        fig.text(1.05, 0.35, 'slope a : {:0.2g} m '.format(a) , fontsize=14, transform=plt.gcf().transFigure)
        fig.text(1.05, 0.3, 'intercept t0 : {:0.2g} m'.format(t0) , fontsize=14, transform=plt.gcf().transFigure)
        fig.text(1.05, 0.25, 'mean residual : {:0.2g} m'.format(mr) , fontsize=14, transform=plt.gcf().transFigure)
        fig.text(1.05, 0.2, '2 standard deviation : {:0.2g} m'.format(sr) , fontsize=14, transform=plt.gcf().transFigure)        

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
        fig.savefig('ths_rapport.png', bbox_inches = 'tight')
    






















