#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 08:56:53 2017

@author: pianarol
"""

import hypy as hp
import scipy as sp
import numpy as np
import math
import matplotlib.pyplot as plt



###function dim ###
def dim(p,t):
    '''THC_DIM - Theis (1941) model with a constant head boundary.

 Syntax: s = thc_dim(p,t)

   p(1) = a  = slope of Jacob Straight Line
   p(2) = t0 = intercept of the first segment of straight line
   p(3) = ti = time of intersection between the 2 straight lines

 Description:
   Computes the drawdown at time t for a constant rate pumping test in 
   a homogeneous confined aquifer bounded by a constant head boundary.

  Reference: Theis, C.V., 1941. The effect of a well on the flow of a
  nearby stream. Transactions of the American Geophysical Union, 22(3):
  734-738.  

 Example:
   s = thc_dim(p,t)

 See also: thc_dmo, thc_gss, thc_rpt'''
    a = p[0]
    
    s = hp.ths.dim([a,p[1]],t) - hp.ths.dim([a,p[2]],t)
    
    return s


###function dls ###
    


def dls(p,t):
    '''THC_DLS - Theis dimensionless drawdown with a constant head boundary

 Syntax: s = thc_dls(p,t)

   p(1) = r1/r2  
   with r1 = radius to the pumping well 
   and  r2 = radius to the image well 

   provides the dimensionless drawdown at reduced time t

 See also: thc_dim, thc_lap
'''
    t = 0.25/t
    
    s = 0.5*sp.special.expn(1,t)-0.5*sp.special.expn(1,t*p[0]**2)
    
    return s
    

###function drw ###


def drw():
    '''THC_DRW - Type curves of the Theis model with a constant head boundary

 Syntax: thc_drw()

 Description:
   Draw a series of type curves of Theis (1935)

 See also: thc_dim, thc_dls
    '''
    
#First figure    
    
    td = np.logspace(-1, 4)
    
    R1 = [1.5,3,10,30]
    R2 = [3,1.5]
    R3 = [10,0]
    R4 = [30,0]
    
    s1 = hp.thc.dls(R1,td)
    s2 = hp.thc.dls(R2,td)
    s3 = hp.thc.dls(R3,td)
    s4 = hp.thc.dls(R4,td)
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

    s = hp.thc.dls(Rd,td)

    xd,yd = hp.ldiff(td,s)
    tdj = hp.jcb.dls(td)
    
    tdrd = []
    
    for i in range(0,len(td)):
        tdrd.append(math.log(Rd[0]))
    
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
    
    
###function frd ###



def frd(rd):
    '''THC_FRD - Function f(rd)

 Syntax: frd = thc_frd(rd)

   rd  = Dimensionless radius
   frd = 2 ln(rd) / ( rd^2 -1 ) rd^(2 rd^2/(1-rd^2))'''
    
    if type(rd) == int:
        if rd == 1:
            rd = 1.00001
        
        rd2 = rd**2
        rd21 = rd2-1
    
        frd = math.log(rd2)/(rd21*math.pow(rd,-2*rd2/rd21))
        
        return frd

    if type(rd) == float:
        if rd == 1:
            rd = 1.00001
        
        rd2 = rd**2
        rd21 = rd2-1
    
        frd = math.log(rd2)/(rd21*math.pow(rd,-2*rd2/rd21))
        
        return frd        

    elif type(rd) == list :
        rd2 = np.power(rd,2)
        rd21 = []
        for i in range(0,len(rd2)):
            rd21.append(rd2[i]-1)
        
        frd = []
        
        rdpow = np.power(rd, -2*rd2/rd21)
        
        rd2log = np.log(rd2)
        frd = np.divide(rd2log,np.multiply(rd21,rdpow))
            
        
        return frd    
    else :
        print('Error, please enter a valid variable (list, int or float)')
    
###function gss ###

def gss(t,s):
    '''THC_GSS - First guess for the parameters of the Theis model with a constant head boundary.

 Syntax:  p = thc_gss(t,s)

   p(1) = a  = slope of Jacob straight line 
   p(2) = t0 = iintercept of the first segment of straight line 
   p(3) = ti = time of intersection between the 2 straight lines

   t    = time
   s    = drawdown

 Description:
   First guess for the parameters of theis solution with a constant head
   boundary.

 See also: thc_dim, thc_dmo, thc_rpt
  '''
    
    td,d = hp.ldiffs(t,s,npoints=10)
 
    sl = np.amax(s)
    sl = np.float(sl)

    du = np.amax(d)
    du = np.float(du)
    i = np.argmax(d)

    
    tu = td[i]
    tu = np.float(tu)

    
    condition = np.greater(t,tu)
    l = np.where(condition == True)
    w = l[0][0]
    
    
    
    
    ta = t[w]

    
    sa = s[w]
    
    
    ts = td[-1]
    ds = d[-1]
    
    #Calculation of the parameters of the model
    
    #Dimensionless distance

    a = sl/du
    

    
    rd = ird(a)

    
    #Slope of Jacob's straight line
    b = sl/2/math.log10(rd)
    
    
    #Origin of jacob straight line
    if rd < 50:
        t0 = 2*tu*math.log(rd)/(0.5625*(math.pow(rd,2)-1))
    else :
        t0 = ta*math.pow(10,-sa/b)
    
    #Origin of the second line
    t1 = math.pow(rd,2)*t0
    
    
    #Time of intersection of the two lines
    ti = math.pow(t1,2)/t0
    
    
    p = []
    
    p.append(b)
    p.append(t0)
    p.append(ti)
    
    return p
    
###function ird ###

def ird(frd1):
    '''THC_IRD - Function inverse of f(rd)

Syntax: rd = thc_ird(frd)

rd  = Dimensionless radius
frd = 2 ln(rd) / ( rd^2 -1 ) rd^(2 rd^2/(1-rd^2))'''    
    
    if type(frd1) == int :
        if frd1 < 2.71828182850322 :
            print('Problem in the inversion of Rd: thc_ird')
            rd = 1.000001
            return rd
        else :
            rd = math.exp(frd1/2)
            
            if rd < 50 :
                y1 = np.arange(1,60.005,0.05)
                y1[0] = 1.000001
                y = []
                for i in range(0,len(y1)):
                    y.append(y1[i])
                
                x = frd(y)
                rd = sp.interpolate.interp1d(x,y)(frd1)
                return rd
            else :
                return rd
            
    if type(frd1) == float :
        if frd1 < 2.71828182850322 :
            print('Problem in the inversion of Rd: thc_ird')
            rd = 1.000001
            return rd
        else :
            rd = math.exp(frd1/2)
            
            
            if rd < 50 :
                y1 = np.arange(1,60.005,0.05)
                y1[0] = 1.000001
                y = []
                for i in range(0,len(y1)):
                    y.append(y1[i])
                
                x = frd(y)
                rd = sp.interpolate.interp1d(x,y)(frd1)
                
                return rd
            else :
                return rd

    else :
        print('Error. Please enter a int or a float')


###function p2x ###

def p2x(p):
    '''%THC_P2X - Conversion of parameters for constant head case

 Syntax: x = thc_p2x(p)

  p(1) = sl  = late time drawdown at the plateau
  p(2) = du  = maximum derivative 
  p(3) = tu  = time of the maximum

  x(1) = a   = slope of Jacob straight line
  x(2) = t0  = intercept of the first Jacob straight line
  x(3) = t1  = intercept of the second Jacob straight line

 Description:


 See also:
'''
    sl = p[0]
    du = p[1]    
    tu = p[2]
    
    sl = np.float(sl)
    du = np.float(du)
    tu = np.float(tu)
    
    
    a = sl/du
    
    rd = ird(a)
    
    x = []
    
    x.append(0.5*sl/math.log(rd))
    x.append(2/0.5626*tu*math.log(rd)/(math.pow(rd,2)-1))
    x.append(math.pow(rd,2)*x[1])
    
    return x
    
###function std ##

def std(l,T,r0):
    '''THC_STD - Computes discharge rate for a well close to a constant head boundary

 Syntax: q = thc_std( l, T, r0 )

   l  = Distance to the hydraulic boundary
   T  = Transmissivity
   r0 = radius of the well

   Nota : the head difference is supposed to be l

 Description:
   The flow field is a half infinite space with a constant head boundary.
   The aquifer is supposed to be confined, homogeneous and the well fully
   penetrating the aquifer.

   The calculation is based on the Dupuit's solution with an image well.
  This equation is also known as the Goodman formula (1965).

 Reference:
   Goodman, R., 1965. Groundwater inflows during tunnel driving.
   Engineering Geology, 2(2): 39-56. 

 Example:
   q=thc_std(100,1e-2,10)

See also: lei_std, thn_std, thc_dim, thc_dmo'''

    q = 2*math.pi*T*l/(math.log(2*l/r0))
    
    return q


###function rpt ###

def rpt(p,t,s,d, name, ttle = 'Interference test', Author = 'My name',  Rapport = 'My Rapport', filetype = 'img'):
    '''%THC_RPT - Report the results of an interpretation with Theis (1941) 
     Syntax: thc_rpt( p, t, s, d, ttle ) 
      p = parameters of the model 
      t = measured time 
      s = measured drawdown 
      d(1) = Q = Pumping rate 
      d(2) = r = Distance between the observation and the pumping well 
      ttle = Title of the figure 
      
      Description: 
          Produces the final figure and results for an interpretation with the 
          Theis model (1941) with a constant head boundary. 
          
          Reference: Theis, C.V., 1941. The effect of a well on the flow of a 
          nearby stream. Transactions of the American Geophysical Union, 22(3): 734-738. 
          
    ssee also: thc_dmo, thc_dim, thc_gss'''

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
    

    tc,sc,mr,sr,rms = hp.script.cmp(p,t,s,'thc')
    
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
        fig.savefig('thc_rapport.pdf', bbox_inches = 'tight')
        
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
        fig.savefig('thc_rapport.png', bbox_inches = 'tight')    

    
    












   
    
    
    
