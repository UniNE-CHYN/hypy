#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 27 17:53:21 2017

@author: pianarol
"""

import hypy as hp
import scipy as sp
import numpy as np
import math
import matplotlib.pyplot as plt

###function lap ###

def lap(x,p):
    '''WAR_LAP - Warren and Root (1963) solution in  Laplace domain

 Syntax:
   war_lap(x,p) provides the dimensionless drawdown at the well

   x(1) = sigma
   x(2) = lamda
   p    = Laplace parameter

 See also: war_dls'''
 
    s = x[0]
    l = x[1]
    
    s = 1/p*sp.special.kv(0,math.sqrt(p+(l*s*p)/(s*p+1)))
    
    return s


###function dls ###
    
def dls(x,td):
    '''WAR_DLS - Dimensionless drawdown of Warren and Root (1963) solution

 Syntax: war_dls(x,t)
   x(1) = sigma 
   x(2) = lamda
   t = time

 Description:
   Calculates the dimensionless drawdown

 References:
   Warren, J. E., and P. J. Root (1963), The behaviour of naturally 
   fractured reservoirs, Society of Petroleum Engineers Journal, 3, 
   245-255.

 See also: war_lap'''
 
    sd = []
    for i in range(0, len(td)):
        sd.append(hp.stefhest('war.lap',x,td[i])) 
        
    return sd



###function dim ###
    
def dim(p,t):
    '''WAR_DIM - Warren and Root (1965) solution

 Syntax: s = war_dim( p, t)

   p(1) = a  = slope of Jacob Straight Line
   p(2) = t0 = intercept with the horizontal axis for 
               the early time asymptote
   p(3) = t1 = intercept with the horizontal axis for 
               the late time asymptote
   p(4) = tm = time of the minimum of the derivative

 Description:
   Calculate the drawdown at time t for confined aquifer with double 
   porosity 

 Example:
   s=war_dim(p,t)

 See also: war_dmo, war_rpt, war_gss'''
 
    a = p[0]
    t0 = p[1]
    t1 = p[2]
    tm = p[3]
    
    td = 0.445268*t/t0
     
    sigma = (t1-t0)/t0
    
    lambada = 2.2458394*t0*math.log(abs(t1/t0))/tm #changed the name, otherwise gives an error
    
    sd = dls([sigma,lambada],td)
    

    
    s = []
    
    for i in range(0,len(sd)):
        s.append(0.868589*a*sd[i])

    return s

###function gss ###
    
def gss(t,s):
    '''WAR_GSS - First guess for the parameters of the Warren and Root solution

 Syntax: p = war_gss(t,s)

   p(1) = a  = slope of Jacob Straight Line
   p(2) = t0 = intercept with the horizontal axis for 
               the early time asymptote
   p(3) = t1 = intercept with the horizontal axis for 
               the late time asymptote
   p(4) = tm = time of the minimum of the derivative

   t    = time
   s    = drawdown

 Description: 
   First guess for the parameters of the Warren and Root solution

   See also: war_dmo, war_dim, war_rpt'''


    td,d = hp.ldiffs(t,s,npoints=40)
    
    dd = np.mean(d[len(d)-4:len(d)])
    
    a = math.log(10)*dd
    
    t0 = t[0]/math.exp(s[-1]/dd)
    
    t1 = t[-1]/math.exp(s[-1]/dd)
    
    
    i = np.argmin(d)
    

    
    
    tm = td[i-2]
    tm = np.float(tm)

    p = []
    
    p.append(a)
    p.append(t0)
    p.append(t1)
    p.append(tm)
    
    return p


###function rpt ###
    
def rpt(p,t,s,d, name, ttle = 'Interference test', Author = 'My name',  Rapport = 'My Rapport', filetype = 'img'):
    '''WAR_RPT - Produces the final figure and results for the Warren and Root model
    
    Syntax: war_rpt( p, t, s, d, ttle ) 
    p(1) = a = slope of Jacob Straight Line 
    p(2) = t0 = intercept with the horizontal axis for the early time asymptote 
    p(3) = t1 = intercept with the horizontal axis for the late time asymptote 
    p(4) = tm = time of the minimum of the derivative 
    
    t = measured time 
     s = measured drawdown 
     d(1) = Q = Pumping rate 
     d(2) = r = Distance to the pumping well 
     ttle = Title of the figure 
     
     Description: 
       Produces the final figure and results for the Warren and Root model
         
       See also: war_dmo, war_dim, war_gss'''
       
   #Rename the parameters for a more intuitive check of the formulas
    Q = d[0]
    r = d[1]
    a = p[0]
    t0 = p[1]
    t1 = p[2]
    tm = p[3]
    
    #Compute the transmissivity, storativity and radius of influence
    Tf=0.1832339*Q/a
    Sf=2.245839*Tf*t0/r**2
    Sm=2.245839*Tf*t1/r**2-Sf;
    
    sigma = (t1-t0)/t0
    lambada = 2.2458394*t0*math.log(t1/t0)/tm
    
    #Calls an internalscript that computes drawdown, derivative and residuals
    #script rpt.cmp    
    

    tc,sc,mr,sr,rms = hp.script.cmp(p,t,s,'war')

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
        fig.text(0.135, -0.1, 'Discharge rate : {:3.2e} m³/s'.format(Q), fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.15, 'Radial distance : {:0.2g} m '.format(r), fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.125, -0.25, 'Hydraulic parameters :', fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.30, 'Transmissivity Tf : {:3.1e} m²/s'.format(Tf), fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.35, 'Storativity Sf : {:3.1e}'.format(Sf), fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.40, 'Storativity Sm : {:3.1e}'.format(Sm), fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.45, 'Interporosity flow lambda : {:0.2e}'.format(lambada), fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.125, -0.55, 'Fitting parameters :' , fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.60, 'slope a : {:0.2g} m '.format(a) , fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.65, 'intercept t0 : {:0.2g} s'.format(t0) , fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.70, 'intercept t1 : {:0.2g} s'.format(t1) , fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.75, 'Minimum deruvatuve tm : {:0.2g} s'.format(tm) , fontsize=14, transform=plt.gcf().transFigure)        



        ax1 = fig.add_subplot(111)
        ax1.set_xlabel('Time in seconds')
        ax1.set_ylabel('Drawdown in meters')
        ax1.set_title(ttle)
        
        ax1.loglog(t, s, c='b',marker = 'o', linestyle = '', label = 'drawdown')
        ax1.loglog(td, sd, c='r',marker = 'x', linestyle = '', label = 'Derivative')
        ax1.loglog(tc, sc, c='g',   label = 'Warren and Root (1965) Model')
        ax1.loglog(tdc, sdc, c='y',  label = 'Model derivative')
        ax1.grid(True)

        ax1.legend()
    
        plt.show()
        fig.savefig('war_rapport.pdf', bbox_inches = 'tight')
        
    if filetype == 'img':
        fig = plt.figure()
        fig.set_size_inches(8, 6)
        fig.text(0.125, 1, Author, fontsize=14, transform=plt.gcf().transFigure) 
        fig.text(0.125, 0.95, Rapport, fontsize=14, transform=plt.gcf().transFigure)  
        
        fig.text(0.125, -0.05, 'Test Data : ', fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.1, 'Discharge rate : {:3.2e} m³/s'.format(Q), fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.15, 'Radial distance : {:0.2g} m '.format(r), fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.125, -0.25, 'Hydraulic parameters :', fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.30, 'Transmissivity Tf : {:3.1e} m²/s'.format(Tf), fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.35, 'Storativity Sf : {:3.1e}'.format(Sf), fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.40, 'Storativity Sm : {:3.1e}'.format(Sm), fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.45, 'Interporosity flow lambda : {:0.2e}'.format(lambada), fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.125, -0.55, 'Fitting parameters :' , fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.60, 'slope a : {:0.2g} m '.format(a) , fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.65, 'intercept t0 : {:0.2g} s'.format(t0) , fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.70, 'intercept t1 : {:0.2g} s'.format(t1) , fontsize=14, transform=plt.gcf().transFigure)
        fig.text(0.135, -0.75, 'Minimum deruvatuve tm : {:0.2g} s'.format(tm) , fontsize=14, transform=plt.gcf().transFigure)
       

        ax1 = fig.add_subplot(111)
        ax1.set_xlabel('Time in seconds')
        ax1.set_ylabel('Drawdown in meters')
        ax1.set_title(ttle)
        
        ax1.loglog(t, s, c='b',marker = 'o', linestyle = '', label = 'drawdown')
        ax1.loglog(td, sd, c='r',marker = 'x', linestyle = '', label = 'Derivative')
        ax1.loglog(tc, sc, c='g',   label = 'Warren and Root (1965) Model')
        ax1.loglog(tdc, sdc, c='y',  label = 'Model derivative')
        ax1.grid(True)

        ax1.legend()
    
        plt.show()
        fig.savefig('war_rapport.png', bbox_inches = 'tight')    
         














