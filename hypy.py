#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 14:48:39 2017

@author: pianarol
"""
import numpy as np
import matplotlib.pyplot as plt
import math
import hypy as hp
from scipy.optimize import leastsq
from scipy import interpolate as spi

### function edit ###  

def edit(file) : 
    '''just show the values of your file'''
    x = np.loadtxt(file, usecols = 0)
    y = np.loadtxt(file, usecols = 1)
    print(x)
    print(y)
    return


### function ldf ###


def ldf(file) :
    '''LDF Load a data file and remove points such that t<=0

 Syntax: [t,s] = ldf( 'fname.dat' )

   fname = filename
   t     = time vector
   s     = drawdown vector

 Description: 
   ldf('filename.dat')is a hytool function designed for loading of data.
   It imports the first and the second column of the file “filename.dat” 
   into the variables t and s (p.e. time and drawdown).

 Example: 
   [t1,s1]=ldf('ths_ds1.dat')

 See also: trial, diagnostic, fit, ths_dmo
'''
    
    t = np.loadtxt(file, usecols = 0)
    #make sure the array contains only float
    t = np.array(t)
    t = np.asfarray(t, float)

#take only the positive values
    condition = np.greater(t,0)
    
    t = np.extract(condition,t)
    
    
    
    s = np.loadtxt(file, usecols = 1)
    
    #make sure the array contains only float
    s = np.array(s)
    s = np.asfarray(s, float)

#take only the positive values
    s = np.extract(condition,s)
    
    return t,s

### function plot ###

def plot(x,y):
    '''create a simple plot with t and s as entry'''
    fig = plt.figure()
            
    ax1 = fig.add_subplot(111)

    ax1.set_title("The Fetter Data set")    #define the title
    ax1.set_xlabel('Time in seconds')       #define the xlabel
    ax1.set_ylabel('Drawdown in meters')    #define the ylabel
                                                 
    ax1.plot(x,y, c='r', label='the data') #create the plot
            
    ax1.legend()                            #add the legend
            
    plt.show()                              #show the plot





##############Differential functions ##########################

###function ldiff ###


def ldiff(t,s):
    '''LDIFF - Approximate logarithmic derivative with centered differences

 Syntax: [xd,yd]=ldiff(x,y)

 See also: ldf, ldiffs, ldiffb, ldiffh'''

#Calculate the difference
    dx = np.diff(t)
    dy = np.diff(s)

#calculate the point

#xd 
    xd = []

    for i in range(0, len(t)-1):
        xd.append( math.sqrt(t[i]*t[i+1]))

#yd
    
    yd = xd*(dy/dx)

    return xd,yd

### function ldiff_plot(t,s)


def ldiff_plot(t,s):
    '''Makes the plot of the ldiff function'''
    xd,yd = ldiff(t,s)
    
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    ax2.set_xlabel('Time')
    ax2.set_ylabel('Drawdown and log derivative')
    
    ax2.loglog(t, s, c='r', label = 'drawdown')
    ax2.scatter(xd, yd, c='b', label = 'derivative')
    ax2.grid(True)

    ax2.legend()
    
    plt.show()


    
### function ldiffs ###


def ldiffs(t,s, npoints = 20):
    '''#LDIFFS - Approximate logarithmic derivative with Spline

 Syntax: [xd,yd] = ldiffs(x,y,[d])
 
   d = optional argument allowing to adjust the number of points 
        used in the Spline

 See also: ldf, ldiff, ldiffb, ldiffh'''

    f = len(t)
#xi and yi
#changing k,s and w affects the graph, make it better or worse

    xi = np.logspace(np.log10(t[0]), np.log10(t[f-1]),  num = npoints, endpoint = True, base = 10.0, dtype = np.float64)

    spl = spi.UnivariateSpline(t,s, k = 5, s = 0.0099)
    yi = spl(xi)


    xd = xi[1:len(xi)-1]
    yd = xd*(yi[2:len(yi)+1]-yi[0:len(yi)-2])/(xi[2:len(xi)+1]-xi[0:len(xi)-2])
    
    return xd,yd

### function ldiffs_plot


def ldiffs_plot(a,b):
    '''Makes the plot of the ldiffs function'''
    xd,yd = ldiffs(a,b)    
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Drawdown and log derivative')
    
    ax1.loglog(a, b, c='r', label = 'drawdown')
    ax1.scatter(xd, yd, c='b', label = 'derivative')
    ax1.grid(True)

    ax1.legend()
    
    plt.show()


###function ldiffb ###


def ldiffb(t,s, d = 2):
    '''#LDIFFB - Approximate logarithmic derivative with Bourdet's formula

 Syntax: [xd,yd]=ldiffb(x,y[,d])

   d = optional argument allowing to adjust the distance between 
       successive points to calculate the derivative.

 See also: ldf, ldiff, ldiffs, ldiffh'''
    ###transform the value of the array X into logarithms 
    logx = []

    for i in range(0, len(t)):
        logx.append(math.log(t[i]))


    dx = np.diff(logx)
    dx1 = dx[0:len(dx)-2*d+1]
    dx2 = dx[2*d-1:len(dx)] 


    dy = np.diff(s)
    dy1 = dy[0:len(dx)-2*d+1]
    dy2 = dy[2*d-1:len(dy)]


#xd and yd

    xd = t[2:len(s)-2]
    yd = (dx2*dy1/dx1+dx1*dy2/dx2)/(dx1+dx2)
    
    return xd,yd
    
###function ldiffd_plot ###


def ldiffb_plot(t,s):
    '''Makes the plot of the ldiffb function'''
    xd,yd = ldiffb(t,s)    

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Drawdown and log derivative')
    
    ax1.loglog(t, s, c='r', label = 'drawdown')
    ax1.scatter(xd, yd, c='b', label = 'derivative')
    ax1.grid(True)

    ax1.legend()
    
    plt.show()
    

###function ldiffh ###



def ldiffh(t,s):
    '''#LDIFFH - Approximate logarithmic derivative with Horne formula

 Syntax: [xd,yd]=ldiffh(t,s)

 See also: ldf, ldiff, ldiffb, ldiffs'''
    #create the table t1,t2,t3 and s1,s2,s3

    endt = len(t)
    ends = len(s)

    t1 = t[0:endt-2]
    t2 = t[1:endt-1]
    t3 = t[2:endt]

    endt1 = len(t1)

    s1 = s[0:ends-2]
    s2 = s[1:ends-1]
    s3 = s[2:ends]


######from the hytool of matlab, to know what is what ##################
#d1 = (log(t2./t1).*s3)./       (log(t3./t2).*log(t3./t1));
#d2 = (log(t3.*t1./t2.^2).*s2)./(log(t3./t2).*log(t2./t1));
#d3 = (log(t3./t2).*s1)./ (log(t2./t1).*log(t3./t1));
#     

#### d1 ####

#log(t2/t1) 
    logt2t1 = []
    for i in range(0, endt1):
        logt2t1.append(math.log(t2[i]/t1[i]))

#log(t2/t1)*s3 : 
    D1_part1 = np.array(logt2t1) * np.array(s3)
 
#log(t3/t2)
    logt3t2 = []
    for i in range(0, endt1):
        logt3t2.append(math.log(t3[i]/t2[i]))

#log(t3/t1)
    logt3t1 = []
    for i in range(0, endt1):
        logt3t1.append(math.log(t3[i]/t2[i]))

#log(t3/t2)*log(t3/t1)
    D1_part2 = np.array(logt3t2)*np.array(logt3t1)
    
    d1 = D1_part1/D1_part2

#### d2 ####

#log(t3*t1/t2²)
    logt3t1t2 = []
    for i in range(0, endt1):
        logt3t1t2.append(math.log(t3[i]*t1[i]/t2[i]**2))

#logt3t1t2 * s2
    D2_part1 = np.array(logt3t1t2) * np.array(s2)

#log(t3/t2)*log(t2/t1)
    D2_part2 = np.array(logt3t2)*np.array(logt2t1)

    d2 = D2_part1 / D2_part2

#### d3 ####

#logt3/t2 * s1
    D3_part1 = np.array(logt3t2) * np.array(s1)

#log(t2/t1)*log(t3/t2)
    D3_part2 = np.array(logt2t1)*np.array(logt3t2)

    d3 = D3_part1 / D3_part2

#xd and yd

    xd = t2
    yd = d1+d2-d3
    
    return xd,yd
#function ldiffh_plot ###

    
def ldiffh_plot(t,s):
    '''Makes the plot of the ldiffh function'''
    xd,yd = ldiffh(t,s)    
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Drawdown and log derivative')
    
    ax1.loglog(t, s, c='r', label = 'drawdown')
    ax1.scatter(xd, yd, c='b', label = 'derivative')
    ax1.grid(True)

    ax1.legend()
    
    plt.show(t,s)


###function diagnostic ###


def diagnostic(a,b, method = 'spline', npoints = 20, step = 2):
    '''DIAGNOSTIC Creates a diagnostic plot of the data
Syntax: diagnostic(t,s, mehtod,npoints,step)

npoints = optional argument allowing to adjust the 
number of points used in the Spline by default
or the distance between points depending on the option below

method = optional argument allowing to select
different methods of computation of the derivative

'spline' for spline resampling
in that case d is the number of points for the spline

'direct' for direct derivation
in that case the value provided in the variable d is not used

'bourdet' for the Bourdet et al. formula
in that case d is the lag distance used to compute the derivative

Description:
This function allows to create rapidly a diagnostic plot
(i.e. a log-log plot of the drawdown as a function of time together with its logarithmic derivative) of the data. 

Example: 
    diagnostic(t,s) 
    diagnostic(t,s,30) 
    diagnostic(t,s,20,'d') - in that case the number 20 is not used
             
See also: trial, ldf, ldiff, ldiffs, fit '''

    if method == 'spline':
        ldiffs_plot(a,b)
    elif method == 'direct' : 
        ldiff_plot(a,b)
    elif method == 'bourdet':
        ldiffb_plot(a,b, d = step)
    else : 
        print('ERROR: diagnostic(t,s,method, number of points)')
        print(' The method selected for log-derivative calculation is unknown')
    
    
##function hyclean ###


def hyclean(t,s):
    '''HYCLEAN - Take only the values that are finite and strictly positive time
                
 Syntax: [tc,sc] = hyclean( t,s )

 Description:
   Take only the values that are finite and strictly positive time 
          
 Example:
   [tc,sc] = hyclean( t,s )

 See also: hyselect, hysampling, hyfilter, hyplot'''
    
    condition = np.logical_and(np.isfinite(s), np.greater(s,0))
    s = np.extract(condition,s)
    t = np.extract(condition, t)
    
    return t,s

###function trial ###
    

def trial(x,t,s, name):
    '''TRIAL Display data and calculated solution together

       Syntax:
           trial(x, t, s, name )          
   
          name = name of the solution
          x    = vector of parameters
          t,s  = data set

       Description:
           The function trial allows to produce a graph that superposes data
           and a model. This can be used to test graphically the quality of a
           fit, or to adjust manually the parameters of a model until a
           satisfactory fit is obtained.

       Example:
           trial(p,t,s,'ths')
           trial([0.1,1e-3],t,s, 'cls')


       See also: ldf, diagnostic, fit, ths_dmo'''
    t,s = hyclean(t,s)
    td,sd = ldiffs(t,s, npoints=40)
    
    tplot = np.logspace(np.log10(t[0]), np.log10(t[len(t)-1]),  endpoint = True, base = 10.0, dtype = np.float64)
    
    if name == 'ths' :
        sc = hp.ths.dim(x,tplot)
    if name == 'Del' :
        sc = hp.Del.dim(x,tplot)
    
    
    tdc,dsc = ldiff(tplot,sc)

    
    if np.mean(sd) < 0 :
        sd = [ -x for x in sd]
        dsc = [ -x for x in dsc]
    condition = np.greater(sd,0)
    td = np.extract(condition,td)
    sd = np.extract(condition,sd)
    
    condition2 = np.greater(dsc,0)
    tdc = np.extract(condition2,tdc)
    dsc = np.extract(condition2,dsc)    
    
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('t')
    ax1.set_ylabel('s')
    ax1.set_title('Log Log diagnostic plot')
    ax1.loglog(t, s, c='r', marker = 'o', linestyle = '')
    ax1.loglog(td,sd, c = 'b', marker = 'x', linestyle = '')
    ax1.loglog(tplot,sc, c = 'g', linestyle = ':')
    ax1.loglog(tdc,dsc, c = 'y', linestyle = '-.')
    
    ax1.grid(True)

    plt.show()
    

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('t')
    ax1.set_ylabel('s')
    ax1.set_title('Semi Log diagnostic plot')
    ax1.semilogx(t, s, c='r', marker = 'o', linestyle = '')
    ax1.semilogx(td,sd, c = 'b', marker = 'x', linestyle = '')
    ax1.semilogx(tplot,sc, c = 'g', linestyle = ':')
    ax1.semilogx(tdc,dsc, c = 'y', linestyle = '-.')
    
    ax1.grid(True)
    
    plt.show()        

###function fit ###



def fit(p0,t,s,name):
    '''FIT - Fit the model parameter of a given model.
#
# Syntax: p = fit(p0, t, s, name)          
#   
#      name   = name of the solution
#      p0     = vector of initial guess for the parameters
#      t,s    = data set
#      
#      p    = vector of the optimum set of parameters
#
# Description: 
#   The function optimizes the value of the parameters of the model so that
#   the model fits the observations. The fit is obtained by an iterative
#   non linear least square procedure. This is why the function requires an
#   initial guess of the parameters, that will then be iterativly modified
#   until a local minimum is obtained.
#   
# Example:
#   p=fit(p0,t,s,'ths')
#
# See also: ldf, diagnostic, trial'''
    
    def residual(p0,t,s,name):
        if name == 'ths':
            sc = hp.ths.dim(p0,t)
        if name == 'Del' :
            sc = hp.Del.dim(p0,t)
            
        scs = []
        
        for i in range(0, len(s)):
            scs.append((sc[i]-s[i])**2)
            
        return scs
    
    def fit2(residual,p0,t,s):
        p, x = leastsq(residual, p0, args = (t,s,name))
    
        return p
    
    p = fit2(residual, p0,t,s)
    
    return p
    
    
    
    
    
    
    
    
    
