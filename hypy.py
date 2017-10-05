#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 14:48:39 2017

@author: pianarol
"""
import numpy as np
import matplotlib.pyplot as plt
import math
import scipy as sp
from scipy import interpolate as spi

##### ldf : use a file and store the data in a vector x,y

def ldf(file) :
    
    x = np.loadtxt(file, usecols = 0)
    #make sure the array contains only float
    x = np.array(x)
    x = np.asfarray(x, float)

#take only the positive values
    inds_t = [i for (i, val) in enumerate(x) if val > 0]
    x = [val for (i, val) in enumerate(x) if i in inds_t]
    
    
    y = np.loadtxt(file, usecols = 1)
    
    #make sure the array contains only float
    y = np.array(y)
    y = np.asfarray(y, float)

#take only the positive values
    inds_t = [i for (i, val) in enumerate(y) if val > 0]
    y = [val for (i, val) in enumerate(y) if i in inds_t]

    return x,y 


###create a simple plot with x and y as entry
def plot(x,y):
    fig = plt.figure()
            
    ax1 = fig.add_subplot(111)

    ax1.set_title("The Fetter Data set")    #define the title
    ax1.set_xlabel('Time in seconds')       #define the xlabel
    ax1.set_ylabel('Drawdown in meters')    #define the ylabel
                                                 
    ax1.plot(x,y, c='r', label='the data') #create the plot
            
    ax1.legend()                            #add the legend
            
    plt.show(x,y)                              #show the plot


#if __name__ == "__main__":              #make the function works, no idea why
#        plot()


##############Differential functions ##########################

#LDIFF - Approximate logarithmic derivative with centered differences
def ldiff(a,b):
######################### For the time t


#Calculate the difference
    dx = np.diff(a)

####################### For the drawdown s

#Calculate the difference
    dy = np.diff(b)

####calculate the point

#xd 
    xd = []

    for i in range(0, len(a)-1):
        xd.append( math.sqrt(a[i]*a[i+1]))


#yd
    yd = xd*(dy/dx)


###plot    
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    ax2.set_xlabel('Time')
    ax2.set_ylabel('Drawdown and log derivative')
    ax2.set_title('Drawdow, Derivative, Location, Northwest')
    ax2.loglog(a, b, c='r', label = 'drawdown')
    ax2.scatter(xd, yd, c='b', label = 'derivative')
    ax2.grid(True)

    ax2.legend()
    
    plt.show(a,b)


    
##ldiffs
def ldiffs(a,b):
######################### For the time t


    f = len(a)
#    a = np.log10(a[f-1])





####################### For the drawdown s

#create w[], contains 1/standard deviation. it's use in the spline function
#give a weight to the points of the graph 
#(see documentation :https://docs.scipy.org/doc/scipy-0.19.1/reference/generated/scipy.interpolate.UnivariateSpline.html )

    ysd = np.std(b)
    ly = len(b)
    w = []

    for i in range(0, ly):
        w.append(1/ysd)


#xi and yi
#changing k,s and w affects the graph, make it better or worse

    xi = np.logspace(np.log10(a[0]), np.log10(a[f-1]),  num = f, endpoint = True, base = 10.0, dtype = np.float64)

    spl = spi.UnivariateSpline(a,b,w, k = 5, s = 0.0099)
    yi = spl(xi)


    xd = xi[1:len(xi)-1]
    yd = xd*(yi[2:len(yi)+1]-yi[0:len(yi)-2])/(xi[2:len(xi)+1]-xi[0:len(xi)-2])
    
    #plot    
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Drawdown and log derivative')
    ax1.set_title('Drawdow, Derivative, Location, Northwest')
    ax1.loglog(a, b, c='r', label = 'drawdown')
    ax1.scatter(xd, yd, c='b', label = 'derivative')
    ax1.grid(True)

    ax1.legend()
    
    plt.show(a,b)


###ldiffh
def ldiffb(a,b, d = 2):
    ###transform the value of the array X into logarithms 
    logx = []

    for i in range(0, len(a)):
        logx.append(math.log(a[i]))


    dx = np.diff(logx)
    dx1 = dx[0:len(dx)-2*d+1]
    dx2 = dx[3:len(dx)] ###########attention au 3, normalement c'est 1*d, mais ne donne pas le résultat attendu

####################### For the drawdown s

    dy = np.diff(b)
    dy1 = dy[0:len(dx)-2*d+1]
    dy2 = dy[3:len(dy)]


#xd and yd

    xd = a[2:len(a)-2]
    yd = (dx2*dy1/dx1+dx1*dy2/dx2)/(dx1+dx2)
    
    #plot    
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Drawdown and log derivative')
    ax1.set_title('Drawdow, Derivative, Location, Northwest')
    ax1.loglog(a, b, c='r', label = 'drawdown')
    ax1.scatter(xd, yd, c='b', label = 'derivative')
    ax1.grid(True)

    ax1.legend()
    
    plt.show(a,b)


###ldiffh
def ldiffh(t,s):
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

    #plot    
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Drawdown and log derivative')
    ax1.set_title('Drawdow, Derivative, Location, Northwest')
    ax1.loglog(t, s, c='r', label = 'drawdown')
    ax1.scatter(xd, yd, c='b', label = 'derivative')
    ax1.grid(True)

    ax1.legend()
    
    plt.show(t,s)

def diagnostic(a,b, m = 's'):
    if m == 's':
        ldiffs(a,b)
    elif m == 'd' : 
        ldiff(a,b)
    elif m == 'b':
        ldiffb(a,b)
    else : 
        print('ERROR: diagnostic(t,s,m)')
        print(' The method selected for log-derivative calculation is unknown')
    
