#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 13:21:11 2017

@author: pianarol
"""

import ldf
import diagnostic
import matplotlib.pyplot as plt

#Define the file you're going to work with

global filename 
filename = '/home/pianarol/Desktop/test.txt'

#We use the ldf function to create to vector, t and s containing the data.
#we first plot them to check graphicaly that they have been correctly loaded. 

def a():
    fig = plt.figure()
        
    ax1 = fig.add_subplot(111)

    ax1.set_title("The Fetter Data set")    #define the title
    ax1.set_xlabel('Time in seconds')       #define the xlabel
    ax1.set_ylabel('Drawdown in meters')    #define the ylabel

    ax1.plot(ldf.t,ldf.s, c='r', label='the data') #create the plot
            
    ax1.legend()                            #add the legend
            
    plt.show()                              #show the plot


if __name__ == "__main__":              #make the function works, no idea why
    a()



##########Diagnostic plot###################

diagnostic.m #Bug with the function. probably declare the variables global! check monday

